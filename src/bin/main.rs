use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::{integration, linear_solvers, advection};
use rust_fluid::boundary::SolidBody;

extern crate image;

fn main() {
    let width = 800;
    let height = 400;

    let mut buffer = vec![0u8; width * height * 3];

    let mut bodies = Vec::new();
    bodies.push(SolidBody::new_box(1.5, 0.7, 0.6, 0.2, std::f64::consts::FRAC_PI_4, 0.0, 0.0, 0.0));
    bodies.push(SolidBody::new_sphere(0.46, 0.5, 0.2, 0.0, 0.0, 0.0, 0.0));

    let mut solver = FluidSolver::new(height, width, 0.005, 1.0 / 400.0, 0.1, bodies)
        .integration(integration::Integration::BogackiShampine)
        .linear_solver(linear_solvers::LinearSolver::ConjugateGradient)
        .advection(advection::Advection::SemiLagrangian);

    // Over 10s
    for iteration in 0..500 {
        for i in 0..4 {
            print!("Step {}: ", 4 * iteration + i);
            solver.add_inflow(0.05, 0.25, 0.2, 0.5, 1.0, 3.0, 0.0);
            solver.update();
        }

        solver.to_image(1.0, &mut buffer);
        image::save_buffer(format!("output/{}.png", iteration), &buffer, width as u32, height as u32, image::RGB(8)).unwrap();
    }
}