"""Inference-only PyTorch backend for the MultiAgeClock R package."""
import gzip
import json
import argparse
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F


def _read_seed(path, layout, device):
    weights = {}
    with gzip.open(path, "rb") as handle:
        for entry in layout:
            shape = tuple(entry["dim"])
            count = int(np.prod(shape))
            raw = handle.read(count * 4)
            if len(raw) != count * 4:
                raise ValueError(f"Incomplete model weight file: {path.name}")
            a = np.frombuffer(raw, dtype="<f4").reshape(shape, order="F")
            if not np.isfinite(a).all():
                raise ValueError(f"Non-finite model weights: {path.name}")
            weights[entry["name"]] = torch.from_numpy(a.copy(order="C")).to(device)
        if handle.read(1):
            raise ValueError(f"Unexpected trailing model data: {path.name}")
    return weights


def _activation(x, name):
    if name == "ELU":
        return F.elu(x)
    if name == "ReLU":
        return F.relu(x)
    if name == "LeakyReLU":
        return F.leaky_relu(x)
    raise ValueError(f"Unsupported activation: {name}")


def _norm(x, w, prefix):
    shape = x.shape
    return F.batch_norm(x.reshape(shape[0], -1),
                        w[prefix + ".running_mean"], w[prefix + ".running_var"],
                        w[prefix + ".weight"], w[prefix + ".bias"],
                        training=False, eps=1e-5).reshape(shape)


def _dense(x, w, prefix, activation, layers, final_linear=False):
    for i in range(layers):
        key = f"{prefix}.layers.{i}"
        x = F.linear(x, w[key + ".weight"], w[key + ".bias"])
        if not (final_linear and i == layers - 1):
            x = _activation(_norm(x, w, f"{prefix}.norms.{i}"), activation)
    return x


def _vectorized(x, w, prefix):
    equation = "bi,hoi->bho" if x.ndim == 2 else "bhi,hoi->bho"
    return torch.einsum(equation, x, w[prefix + ".weight"]) + w[prefix + ".bias"].unsqueeze(0)


def _forward(x, w, activation):
    shared = _dense(x, w, "shared", activation, 3)
    bypass = x
    for i in range(1, 4):
        bypass = _vectorized(bypass, w, f"disease_tower.bypass_{i}")
        bypass = _activation(_norm(bypass, w, f"disease_tower.norm_{i}"), activation)
    shared_heads = shared.unsqueeze(1).expand(-1, bypass.shape[1], -1)
    disease = _vectorized(torch.cat((shared_heads, bypass), dim=2), w, "disease_tower.predictor_1")
    disease = _activation(_norm(disease, w, "disease_tower.predictor_norm"), activation)
    disease = _vectorized(disease, w, "disease_tower.predictor_2").squeeze(-1)
    mortality = _dense(x, w, "mortality_head.bypass", activation, 3)
    mortality = _dense(torch.cat((shared, mortality), dim=1), w,
                       "mortality_head.predictor", activation, 2, final_linear=True)
    return torch.cat((disease, mortality), dim=1)


def predict_scores(values, model_json, model_dir, batch_size=256, device="cpu"):
    """Return five-seed mean raw scores, with no fitting or missing-value imputation."""
    with open(model_json, encoding="utf-8") as handle:
        meta = json.load(handle)
    values = np.asarray(values, dtype=np.float32)
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] != len(meta["features"]):
        raise ValueError("Input matrix does not match the model feature dimensions.")
    if not np.isfinite(values).all():
        raise ValueError("Complete finite inputs are required.")
    batch_size = int(batch_size)
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")
    device = torch.device(device)
    x = ((values - np.asarray(meta["means"], dtype=np.float32)) /
         np.asarray(meta["sds"], dtype=np.float32)).astype(np.float32)
    scores = np.zeros((len(x), len(meta["endpoints"])), dtype=np.float64)
    previous_threads = torch.get_num_threads()
    try:
        if device.type == "cpu":
            torch.set_num_threads(min(previous_threads, 8))
        with torch.inference_mode():
            for seed in meta["seed_files"]:
                w = _read_seed(Path(model_dir) / meta["model_id"] / seed, meta["tensor_layout"], device)
                for start in range(0, len(x), batch_size):
                    batch = torch.from_numpy(np.ascontiguousarray(x[start:start + batch_size])).to(device)
                    scores[start:start + batch_size] += _forward(batch, w, meta["activation"]).cpu().numpy()
                del w
    finally:
        if device.type == "cpu":
            torch.set_num_threads(previous_threads)
    return scores / len(meta["seed_files"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply frozen MultiAgeClock inference weights.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--weights", required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--columns", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    values = np.fromfile(args.input, dtype="<f8").reshape((args.rows, args.columns), order="F")
    scores = predict_scores(values, args.metadata, args.weights, args.batch_size, args.device)
    scores.astype("<f8").tofile(args.output)
