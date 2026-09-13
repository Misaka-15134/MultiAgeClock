# Frozen models and input contract

Each neural model uses a shared representation, 36 disease score heads and one mortality score head. The arithmetic mean of five frozen seed-specific scores is converted to biological age using frozen coefficients. The inference implementation evaluates linear layers, batch normalization in evaluation mode and the saved activation function. Dropout is inactive at inference.

For a single modality, the necessary score coefficients already combine the frozen endpoint coefficient and aggregation weight. Integrated scores are formed from the three full-modality scores using the saved endpoint-specific stacking coefficients.

The age conversion is

\[
\begin{aligned}
\eta &= \sum_j c_j s_j,\\
BA &= CA + a(\eta-\bar\eta),\\
\mathrm{raw\_gap} &= BA-CA,\\
BAA &= BA-f_{\mathrm{reference}}(CA).
\end{aligned}
\]

Here the coefficients, reference mean, age scale and natural cubic spline are frozen. R evaluates the reference spline from its knot locations and fitted values. The package does not fit any quantity to a user's table.

The public model metadata contain only inference parameters, tensor dimensions, input dictionaries and file names. Each compressed weight file stores float32 tensors in little-endian, column-major order matching the declared tensor layout. Model loading checks that every declared tensor is present and finite. Optimizer states, training records, training imputation donors, checkpoint provenance and participant identifiers are not part of these assets.

The public reference examples are generated artificial inputs, with expected results calculated by the original frozen PyTorch CPU implementation. Small floating-point differences between R double-precision matrix operations and PyTorch float32 are expected. Numerical tests use an absolute tolerance of 0.001 years for BA, raw gap and BAA on these examples.

All inputs must be complete. Full clinical and integrated feature names follow the current frozen clinical binding dictionary. The K17 panel is the separately frozen common clinical panel used for external validation. Earlier candidate panels are not substituted for it.

Input unit metadata are available through `model_features()`. Extra columns are ignored, column order does not affect the calculation, and ambiguous mappings are rejected. There is no automatic unit inference, imputation, panel substitution or target-cohort recalibration.
