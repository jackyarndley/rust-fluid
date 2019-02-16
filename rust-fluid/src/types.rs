use crate::field::Field;

pub type IntegrableFunction = Fn(f64, f64) -> f64;
pub type Integration = fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64;
pub type Interpolation = fn(f64, f64, &Field) -> f64;
pub type Advection = fn(&mut Field, &Field, &Field, f64, f64, &Fn(f64, f64, &Field) -> f64, &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64);
pub type LinearSolver = fn(&mut Field, &Field, f64, f64, f64, usize);