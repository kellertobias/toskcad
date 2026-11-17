//! Developer probe of original supplementary STEP curve agreement.
use multicad_rust_kernel::step_read::brep::validate_pcurves;
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let inspect_refusal = args.iter().any(|a| a == "--retained-geometry-refusal");
    let mut failed = false;
    for path in args.iter().filter(|a| *a != "--retained-geometry-refusal") {
        let input = std::fs::read_to_string(path).expect("read STEP input");
        match validate_pcurves(&input, 1e-7) {
            Ok(report) => println!("{path}: {report:?}"),
            Err(error) => {
                println!("{path}: {error:?}");
                failed = true;
            }
        }
        if inspect_refusal {
            println!(
                "{path}: retained_geometry_refusal={:?} (diagnostic only; no admission)",
                multicad_rust_kernel::step_read::brep::retained_geometry_refusal(
                    &input, 1e-7, None
                )
            );
        }
    }
    if failed {
        std::process::exit(1)
    }
}
