use crate::util::Field;
use crate::util::Sparse;

pub type IntegrableFunction = Fn(f64, f64) -> f64;
pub type Integration = fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64;
pub type Interpolation = fn(f64, f64, &Field) -> f64;
pub type Advection = fn(&mut Field, &Field, &Field, &Field, f64, f64, &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64);
pub type LinearSolver = fn(&mut Vec<f64>, &mut Vec<f64>, &mut Vec<f64>, &mut Vec<f64>, &mut Vec<f64>, &mut Sparse, f64, f64, f64, usize, usize, usize);