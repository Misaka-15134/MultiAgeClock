# MultiAgeClock

Calculate frozen multimodal biological ages from complete tabular measurements in R.

| Model ID | Biomarker inputs | Additional input |
|---|---:|---|
| `clinical_k17` | 17 clinical measurements | Age in years |
| `clinical` | 47 clinical measurements | Age in years |
| `olink` | 2,920 Olink proteins | Age in years |
| `nmr` | 168 NMR measurements | Age in years |
| `integrated` | All three full panels, 3,135 measurements | Age in years |

The K17 model is the frozen clinical panel used for external validation. Each neural model averages five frozen seed models. Inference runs locally. The default R backend does not require Python; an optional PyTorch backend supports vectorized CPU and GPU inference in a separate local Python process.

```r
install.packages("remotes")
remotes::install_github("Misaka-15134/MultiAgeClock@v0.1.0")
library(MultiAgeClock)

list_models()
model_features("clinical_k17")
download_models("clinical_k17")  # One-time download, stored in a versioned cache

model <- load_model("clinical_k17")
result <- predict_age(example_data(), model, id_col = "sample_id")
result

# Score a real file; input data remain on your computer.
score_file("measurements.csv", "biological_ages.csv",
           model = model, id_col = "sample_id")
```

`ba` is biological age in years. `raw_gap` is BA minus chronological age. `baa` is the residual after subtracting the frozen reference-age spline, also in years. No model or calibration parameter is fitted to the input table. Values are not clipped to the training age range; `age_outside_reference` flags extrapolation.

To use an existing Python environment with NumPy and PyTorch, pass `backend = "pytorch", python = "/path/to/python"` to `predict_age()` or `score_file()`. Set `device = "cuda"` for a configured CUDA environment. Both backends use the same frozen parameters.

Use the exact units returned by `model_features()`. Inputs must be unstandardized measurements on the documented assay scale. Olink inputs are NPX values already expressed on their log2 scale. NMR names are the frozen UK Biobank field identifiers; their Nightingale names and units are included in the dictionary. Missing columns, missing values and non-finite values produce an error. No training-subject imputation database is distributed.

CSV and TSV are supported directly. For Excel files, install `readxl` for reading and `writexl` for writing. The package preserves row order and text sample identifiers and refuses to overwrite existing outputs.

[中文教程](docs/tutorial_zh.md) · [Model and input notes](docs/models.md) · [Code availability](docs/code_availability.md)

This repository provides frozen-model inference, necessary preprocessing and calibration parameters, synthetic examples and tests. Training, hyperparameter search, feature selection, research analyses and participant data are outside this distribution. The models are intended for research; individual results do not establish clinical diagnostic thresholds.

Model files are distributed separately in [versioned releases](https://github.com/Misaka-15134/MultiAgeClock/releases). Installation from the source repository alone does not include the weight files. See [LICENSE](LICENSE) for usage terms.
