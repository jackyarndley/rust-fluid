use util::Sparse;

mod gauss_siedel;
mod conjugate_gradient;

pub use self::gauss_siedel::*;
pub use self::conjugate_gradient::*;

pub fn empty(_: &mut Vec<f64>, _: &mut Vec<f64>, _: &mut Vec<f64>, _: &mut Vec<f64>, _: &mut Vec<f64>, _: &mut Sparse, _: f64, _: f64, _: f64, _: usize, _: usize, _:usize) {

}

