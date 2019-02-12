use rust_fluid::fluid_solver::FluidSolver;


fn main() {
    let mut buffer = vec![0u8; 128 * 128 * 4];

    let mut solver = FluidSolver::new(128, 128, 0.1);

    let mut time: f64 = 0.0;
    let mut iterations: usize = 0;

    while time < 1.0 {
        solver.add_source2(20, 20, 40, 40, 2.0, 2.0, 1.0);

        solver.update(0.005);

        solver.to_image(&mut buffer);

        lodepng::encode32_file(format!("output/{}.png", iterations), &buffer, 128, 128).unwrap();
        iterations += 1;
        time += 0.005;
    }
}