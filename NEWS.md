# ssmodels 2.0.0

## Major updates
- Complete overhaul of the package, improving organization, readability, and performance of all functions.
- Rewritten log-likelihood and gradient functions (`loglik_*` and `gradlik_*`) for enhanced numerical stability and clarity.
- Fixed discrepancies where analytical gradients did not match numerical gradients.
- Comprehensive documentation updates for all functions, ensuring better understanding and usage.
- Added two new helper functions:
  - `postprocess_theta()`: streamlines parameter transformations for clear interpretation and improved consistency across models.
  - `extract_model_components()`: extracts `model.frame`, `model.matrix`, and `model.response` objects in a robust and reusable way.
- All functions now follow consistent coding style and best practices.
- Significant performance improvements, making the package lighter and more efficient.

## Bug fixes
- Fixed issues with incorrect gradient calculations for `sigma` and `rho` parameters.
- Corrected numerical errors in several model functions.

## Other improvements
- Updated vignette and examples to reflect the new structure and improvements.
- Switched pkgdown site to Bootstrap 5 for improved readability and responsiveness.

# ssmodels 1.0.1

## Minor updates
- Improved documentation and examples.
- Added unit tests to ensure stability of `HeckmanCL()` and other core functions.

# ssmodels 1.0.0

## Initial release
- Initial implementation of the classic Heckman model (`HeckmanCL()`) and foundational sample selection models.
- Basic infrastructure for selection bias correction in econometric models.
