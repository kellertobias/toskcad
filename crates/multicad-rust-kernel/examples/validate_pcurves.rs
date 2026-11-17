//! Developer probe of original supplementary STEP curve agreement.
use multicad_rust_kernel::step_read::brep::validate_pcurves;
fn main() {
    let mut failed = false;
    for path in std::env::args().skip(1) {
        let input = std::fs::read_to_string(&path).expect("read STEP input");
        match validate_pcurves(&input, 1e-7) {
            Ok(report) => println!("{path}: {report:?}"),
            Err(error) => {
                println!("{path}: {error:?}");
                failed = true;
            }
        }
    }
    if failed {
        std::process::exit(1)
    }
}
