//! Emit measured Quasar costs as JSON, one object per configuration, for
//! `quasar_cost_model.py` to diff its predictions against.
//!
//! Usage: `quasar_estimate <nv:log_rows:index>[,...]`
//!   e.g. `cargo run -p pcs-profile --release --bin quasar_estimate -- 14:6:2,16:6:2,16:8:2,16:6:4`
//!
//! Column queries come from the proven QA distance over ed25519 (clamped to the
//! codeword for small instances). Build with the default `count-ops` feature for
//! the arithmetic columns; turn it off for meaningful wall-clock.

use pcs_profile::quasar::QuasarProfile;
use pcs_profile::{PcsMetrics, ProfileablePcs};

fn main() {
    let spec = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "12:6:2,14:6:2,16:6:2".to_string());

    println!("[");
    let configs: Vec<&str> = spec.split(',').collect();
    for (i, cfg) in configs.iter().enumerate() {
        let parts: Vec<usize> = cfg
            .split(':')
            .map(|p| p.trim().parse().expect("nv:log_rows:index"))
            .collect();
        assert_eq!(parts.len(), 3, "each config is nv:log_rows:index");
        let (num_vars, log_rows, index) = (parts[0], parts[1], parts[2]);

        let profile = QuasarProfile {
            log_rows,
            index,
            ..QuasarProfile::default()
        };
        let q = profile
            .security_queries(num_vars)
            .expect("ed25519 admits a QA distance");
        let m: PcsMetrics = profile.run(num_vars, q, 1);

        let log_cols = num_vars - log_rows;
        let bf_queries = profile.bf_queries.min(profile.bf_inv_rate << log_cols);
        print!(
            "  {{\"num_vars\": {}, \"log_rows\": {}, \"index\": {}, \"bf_inv_rate\": {}, \
             \"bf_queries\": {}, \"q\": {}, \"sec_bits\": {:.2}, \"proof_bytes\": {}, \
             \"hash_leaf\": {}, \"hash_compress\": {}, \"hash_sponge\": {}, \
             \"mul\": {}, \"sq\": {}, \"inv\": {}, \"add\": {}, \
             \"commit_us\": {}, \"open_us\": {}, \"verify_us\": {}}}{}",
            num_vars,
            log_rows,
            index,
            profile.bf_inv_rate,
            bf_queries,
            m.q,
            m.sec_bits,
            m.proof_bytes,
            m.hash_leaf,
            m.hash_compress,
            m.hash_sponge,
            m.mul,
            m.sq,
            m.inv,
            m.add,
            m.commit_us,
            m.open_us,
            m.verify_us,
            if i + 1 == configs.len() { "" } else { "," }
        );
        println!();
    }
    println!("]");
}
