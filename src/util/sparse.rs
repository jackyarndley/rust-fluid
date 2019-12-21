pub struct Sparse {
    pub diagonals: Vec<f64>,
    pub plus_x:    Vec<f64>,
    pub plus_y:    Vec<f64>,
}

impl Sparse {
    pub fn new(size: usize) -> Sparse {
        Sparse {
            diagonals: vec![0.0; size],
            plus_x:    vec![0.0; size],
            plus_y:    vec![0.0; size]
        }
    }
}