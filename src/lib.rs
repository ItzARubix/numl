//! # numl
//! Implementations of various numerical algorithms with an emphasis on accuracy. 

use thiserror::Error;

/// Enum of errors that can be returned by numl functions.
#[derive(Error, Debug)]
pub enum NumlError {
    /// Error for when a typical value of zero is passed into a function that requires a
    /// nonzero typical value.
    ///
    /// Functions sometimes require a typical value to be passed in in order to ensure that the
    /// underlying algorithm doesn't fail for small input values. As a result, most algorithms
    /// cannot guarantee accuracy if the typical value is set to zero.
    ///
    /// It is recommended that you also avoid passing typical values with very small magnitude into
    /// these functions, though doing so will not return an error unless the value is exactly zero.
    #[error("Typ must be positive")]
    TypError,

    /// Error returned by functions which require a derivative to be zero. 
    /// This is usually returned when a function would like to divide by a calculated derivative,
    /// in which case a zero derivative would cause an erroneous divide-by-zero.
    #[error("Derivative calculated to zero, but needs to be nonzero")]
    DerivativeZeroError,
}

/// Numerically calculates the derivative of the given function at the specified point.
///
/// Inputs:
/// - f: fn(f64) -> f64
/// - x: f64
/// - typ: f64
///
/// f() is the function whose derivative is being computed, x is the point at which that
/// derivative is computed, and typ is the typical size of x (in the event that the passed value of
/// x is very different from the usual size that x takes on). Please see the documentation of
/// NumlError::TypError for more information on the typical value parameter.
///
/// f() is expected to be a pure function, and this algorithm will do two evaluations of the
/// function in order to determine the derivative.
pub fn derivative(f: fn(f64) -> f64, x: f64, typ: f64) -> Result<f64, NumlError> {

    if typ==0.0 {
        return Err(NumlError::TypError); 
    }

    let h:f64 = if x.abs() > typ.abs() { 
        f64::cbrt(f64::EPSILON)*x
    } else {
        f64::cbrt(f64::EPSILON)*typ
    };

    Ok((f(x+h) - f(x-h))/(2.0*h))
}

/// Performs one iteration of a quasi-Newton's method and returns the result.
///
/// Inputs:
/// - f: fn(f64) -> f64
/// - x: f64
/// - typ: f64
///
/// f() is the function whose root is being computed, x is the current guess of the root, 
/// and typ is the typical size of x (in the event that the passed value of
/// x is very different from the usual size that x takes on). Please see the documentation of
/// NumlError::TypError for more information on the typical value parameter.
///
/// f() is expected to be a pure function, and this algorithm will do three evaluations of the
/// function in order to determine the derivative.
///
/// If the derivative of f() at the specified point is evaluated to be exactly a floating point
/// zero, a NumlError::DerivativeZeroError will be returned. However, no error will be returned if
/// the derivative evalutes to a number very close to zero, which may cause issues. It is thus
/// recommended to check if your function has a derivative zero near the input value if you are
/// getting unexplainable behavior.
pub fn nqn(f: fn(f64) -> f64, x: f64, typ: f64) -> Result<f64, NumlError> { 
    let computed_derivative:f64 = match derivative(f, x, typ) {
        Ok(0.0) => return Err(NumlError::DerivativeZeroError),
        Ok(g) => g,
        Err(g) => return Err(g)
    };
    Ok((x) - (f(x))/(computed_derivative))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_cubic(x: f64) -> f64 {
        (x*x*x) + (2.0*x*x) - 0.4
    }

    #[test]
    fn test_derivative() {
        let result = derivative(sample_cubic, 1.0, 0.5).unwrap();
        assert!(result > 6.9 && result < 7.1);
    }

    #[test]
    fn test_derivative_neg() {
        let result = derivative(sample_cubic, -1.0, 0.5).unwrap();
        assert!(result < -0.9 && result > -1.1);
    }

    #[test]
    fn test_derivative_negtyp() {
        let result = derivative(sample_cubic, 1.0, -0.5).unwrap();
        assert!(result > 6.9 && result < 7.1);
    }

    #[test]
    fn test_derivative_neg_negtyp() {
        let result = derivative(sample_cubic, -1.0, -0.5).unwrap();
        assert!(result < -0.9 && result > -1.1);
    }

    #[test]
    fn test_nqn() {
        let mut guess = 1.0;
        for _i in 1..12 {
            guess = nqn(sample_cubic, guess, 0.5).unwrap();
        }
        assert!(guess > 0.4 && guess < 0.41);
    }
    
}
