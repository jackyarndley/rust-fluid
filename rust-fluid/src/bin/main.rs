use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::advection;
use rust_fluid::interpolation;
use rust_fluid::integration;
use rust_fluid::linear_solvers;


fn main() {
    let mut buffer = vec![0u8; 256 * 512 * 4];

    let mut solver = FluidSolver::new(256, 512, 0.01, 1.0 / 8.0, 1.0)
        .interpolation(interpolation::bicubic_interpolate)
        .integration(integration::bogacki_shampine)
        .linear_solver(linear_solvers::gauss_siedel)
        .advection(advection::semi_lagrangian);

    solver.y_velocity_src.add_inflow(80, 110, 112, 142, 75.0);
    solver.density_src.add_inflow(80, 110, 112, 142, 2.0);

    for iteration in 0..400 {
        solver.solve();

        solver.y_velocity_src.add_inflow(80, 110, 112, 142, 75.0);
        solver.density_src.add_inflow(80, 110, 112, 142, 2.0);

        solver.to_image(&mut buffer);
        lodepng::encode32_file(format!("output/{}.png", iteration), &buffer, 512, 256).unwrap();
    }
}
