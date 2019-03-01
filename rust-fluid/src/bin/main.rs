use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::advection;
use rust_fluid::interpolation;
use rust_fluid::integration;
use rust_fluid::linear_solvers;


fn main() {
    let mut buffer = vec![0u8; 128 * 128 * 4];

    let mut solver = FluidSolver::new(128, 128, 0.005, 1.0 / 128.0, 0.1)
        .interpolation(interpolation::bicubic_interpolate)
        .integration(integration::bogacki_shampine)
        .linear_solver(linear_solvers::gauss_siedel)
        .advection(advection::semi_lagrangian);

    solver.y_velocity_src.add_inflow(58, 26, 77, 30, 1.0);
    solver.density_src.add_inflow(58, 26, 77, 30, 3.0);

    for iteration in 0..1000 {
        for i in 0..4 {
            print!("Global iteration {}. ", 4 * iteration + i);
            solver.solve();

            solver.y_velocity_src.add_inflow(58, 26, 77, 30, 1.0);
            solver.density_src.add_inflow(58, 26, 77, 30, 3.0);
        }


        solver.to_image(&mut buffer);
        lodepng::encode32_file(format!("output/{}.png", iteration), &buffer, 128, 128).unwrap();
    }
}
