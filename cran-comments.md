## Test results

- I have run `devtools::check()` on Windows, macOS, and Linux.  
- I have also run `devtools::check_win_devel()` on the win-builder.  
- I have run `devtools::spell_check()` and all words have been reviewed.  
- I have run `rhub::check()` on multiple platforms.

## Changes in this version (2.0.0)

- Improved the likelihood and gradient functions in all model functions for better speed and stability.  
- Documented all functions comprehensively, following best practices.  
- Added two new helper functions: `postprocess_theta()` and `extract_model_components()`, to promote modular and clean coding.  
- Fixed analytical gradients to ensure they match the numerical gradients in all cases.  
- The package is now lighter, faster, and more functional overall.

## Additional information

- There are no notes, warnings, or errors from the CRAN checks.  
- All issues related to package website (`pkgdown`) formatting have been resolved.  
- The package website is available at: https://fsbmat-ufv.github.io/ssmodels/

Thank you for reviewing my submission!  
Fernando de Souza Bastos
