//! Sweep PCS-only cost (BabyBear hash) over num_vars for Basefold and the three real Brakedown
//! presets (small-proof / balanced / fast-prover), each with a distance-derived query count.
//!
//! Usage: `sweep [basefold|brakedown|quasar|both|all] [num_vars_csv] [prov|conj|both]`
//!   e.g. `cargo run -p pcs-profile --release --bin sweep -- both 8,10,12,14,16 both`
//!
//! Two cost centers per row: HASH perms (leaf/compress/sponge) and ARITH (mul/sq/inv/add), + proof
//! bytes. NOTE: Brakedown's real query count (thousands) exceeds the codeword for small num_vars, so
//! `q` clamps to the codeword length there (fully opened) — the δ/3 vs δ/2 gap only shows once the
//! codeword is large enough (num_vars ≳ 22 for the small-proof preset). Build with `--features
//! count-ops` (default) for the arithmetic columns.

use p3_brakedown::params;
use pcs_profile::basefold::BasefoldProfile;
use pcs_profile::brakedown::BrakedownProfile;
use pcs_profile::quasar::QuasarProfile;
use pcs_profile::{PcsMetrics, ProfileablePcs};

fn main() {
    let which = std::env::args().nth(1).unwrap_or_else(|| "both".to_string());
    let nvs: Vec<usize> = std::env::args()
        .nth(2)
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![8, 10, 12, 14, 16]);
    let regime = std::env::args().nth(3).unwrap_or_else(|| "both".to_string());

    let bf = BasefoldProfile::default();
    let regimes: Vec<bool> = match regime.as_str() {
        "prov" => vec![false],
        "conj" => vec![true],
        _ => vec![false, true],
    };
    let mut bds: Vec<BrakedownProfile> = Vec::new();
    for &conj in &regimes {
        for p in params::all() {
            bds.push(BrakedownProfile::new(p, conj));
        }
    }

    println!(
        "# PCS-only verify cost (hash = BabyBear Poseidon2-16, λ=100). Basefold q=288 (list); \
         Brakedown presets: real t=⌈λ/−log₂(1−e)⌉, reg d/3=provable e=δ/3, d/2=conjectured e=δ/2."
    );
    println!("# (Brakedown q clamps to codeword for small num_vars — real counts show at num_vars ≳ 22.)");
    println!(
        "# Quasar: q = committed columns from the QA distance bound (reg d/3), \
         plus 288 BaseFold reps for the code-switching openings."
    );
    println!("{}", PcsMetrics::header());
    let quasars: Vec<QuasarProfile> = if which == "quasar" || which == "all" {
        // Both rates, at the paper's own benchmark shape (N1 = 2^6).
        vec![
            QuasarProfile::default(),
            QuasarProfile {
                index: 4,
                ..QuasarProfile::default()
            },
        ]
    } else {
        Vec::new()
    };

    for &nv in &nvs {
        if which == "basefold" || which == "both" || which == "all" {
            println!("{}", bf.run(nv, 288, 1).row());
        }
        if which == "brakedown" || which == "both" || which == "all" {
            for bd in &bds {
                println!("{}", bd.run(nv, 0, 1).row());
            }
        }
        for qz in &quasars {
            // Skip shapes with fewer than two columns.
            if nv <= qz.log_rows {
                continue;
            }
            let q = qz.security_queries(nv).expect("ed25519 admits a QA distance");
            println!("{}", qz.run(nv, q, 1).row());
        }
    }
}
