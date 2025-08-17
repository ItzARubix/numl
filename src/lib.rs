//! # numl
//! Implementations of various numerical algorithms with an emphasis on accuracy. 


use thiserror::Error;

/// Enum of errors that can be returned by numl functions.
/// See the #[error(...)] messages for information on particular errors.
#[derive(Error, Debug)]
pub enum NumlError {
    #[error("Typ must be positive")]
    TypError,
    #[error("Derivative calculated to zero, but needs to be nonzero")]
    DerivativeZeroError,
}

/// Numerically calculates the derivative of the given function at the specified point.
pub fn derivative(f: fn(f64) -> f64, x: f64, typ: f64) -> Result<f64, NumlError> {

    if !(typ>0.0) {
        return Err(NumlError::TypError); 
    }

    let h:f64 = if x.abs() > typ { 
        f64::cbrt(f64::EPSILON)*x
    } else {
        f64::cbrt(f64::EPSILON)*typ
    };

    Ok((f(x+h) - f(x-h))/(2.0*h))
}

/// Performs one iteration of a quasi-Newton's method and returns the result.
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
    fn test_nqn() {
        let mut guess = 1.0;
        for _i in 1..12 {
            guess = nqn(sample_cubic, guess, 0.5).unwrap();
        }
        assert!(guess > 0.4 && guess < 0.41);
    }
    
}
