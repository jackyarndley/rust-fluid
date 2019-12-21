

pub struct ScalarGrid2 {
    pub data: Vec<f64>,
    pub rows: usize,
    pub columns: usize,
    pub len: usize,
    pub center: (f64, f64)
}

impl ScalarGrid2 {
    pub fn new(rows: usize, columns: usize) -> ScalarGrid2 {
        ScalarGrid2 {
            data: vec![0.0; rows * columns],
            rows,
            columns,
            len: rows * columns,
            center: (0.5, 0.5)
        }
    }
}
