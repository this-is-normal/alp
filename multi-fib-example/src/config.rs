//! The multi-stark backend wiring: field, hash, Merkle MMCS, and the WHIR
//! multilinear PCS. This is the module that will change when we swap WHIR for
//! the trivial PCS (M1) and later Brakedown (M3) — the AIR and `main` stay put.

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_brakedown::{
    BrakedownCode as RealBrakedownCode, BrakedownMlPcs, LinearCode, TensorPcs,
};
use p3_challenger::DuplexChallenger;
use p3_dft::Radix2DFTSmallBatch;
use p3_field::Field;
use p3_field::extension::BinomialExtensionField;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multi_stark::config::MultiStarkConfig;
use p3_sumcheck::layout::{Layout, PrefixProver, Table, Witness};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_trivial_mpcs::TrivialMultilinearPcs;
use p3_util::log2_ceil_usize;
use p3_whir::{
    DomainSeparator, FoldingFactor, ProtocolParameters, SecurityAssumption, WhirConfig, WhirProver,
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

// ---- field / hash / commitment stack (BabyBear + Poseidon2 + Merkle + WHIR) ----
pub type F = BabyBear;
pub type EF = BinomialExtensionField<F, 4>;
type Perm = Poseidon2BabyBear<16>;

type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
pub type MyChallenger = DuplexChallenger<F, Perm, 16, 8>;

type PackedF = <F as Field>::Packing;
type MyMmcs = MerkleTreeMmcs<PackedF, PackedF, MyHash, MyCompress, 2, 8>;

type MyDft = Radix2DFTSmallBatch<F>;
type L = PrefixProver<F, EF>;
type MyPcs = WhirProver<EF, F, MyDft, MyMmcs, MyChallenger, L>;

/// First-round folding factor; also the per-table padding floor.
const FOLDING: usize = 2;

/// A WHIR-backed multilinear AIR configuration over BabyBear.
pub struct MultiFibConfig {
    pcs: MyPcs,
}

impl MultiStarkConfig for MultiFibConfig {
    type Val = F;
    type Challenge = EF;
    type Challenger = MyChallenger;
    type Pcs = MyPcs;

    fn pcs(&self) -> &MyPcs {
        &self.pcs
    }

    fn min_num_variables(&self) -> usize {
        FOLDING
    }

    fn build_witness(&self, tables: Vec<Table<F>>) -> Witness<F> {
        L::new_witness(tables, FOLDING)
    }

    fn committed_table<'a>(
        &self,
        prover_data: &'a p3_whir::WhirProverData<F, EF, MyMmcs, L>,
        table_index: usize,
    ) -> &'a Table<F> {
        prover_data.table(table_index)
    }
}

/// Fixed permutation so the prover and verifier transcripts match exactly.
fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(0xD15EA5E);
    Perm::new_from_rng_128(&mut rng)
}

/// Per-round log-inverse rates for a stacked polynomial.
fn default_round_log_inv_rates(num_variables: usize, folding_factor: &FoldingFactor) -> Vec<usize> {
    let folding_schedule = folding_factor
        .compute_folding_schedule(num_variables)
        .expect("valid folding schedule");
    let num_rounds = folding_schedule.len() - 1;
    let mut rates = Vec::with_capacity(num_rounds);
    let mut rate = 1;
    for &folding in folding_schedule.iter().take(num_rounds) {
        rate += folding - 1;
        rates.push(rate);
    }
    rates
}

/// Build a configuration sized for a single trace of `2^log_height` rows and `width` columns.
///
/// The committed polynomial stacks every trace column, so WHIR is sized for the
/// stacked variable count `log_height + ceil(log2(width))`.
pub fn config_for(log_height: usize, width: usize) -> MultiFibConfig {
    let stacked_num_variables = log_height + log2_ceil_usize(width);
    let folding_factor = FoldingFactor::Constant(FOLDING);

    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let params = ProtocolParameters {
        security_level: 32,
        pow_bits: 0,
        round_log_inv_rates: default_round_log_inv_rates(stacked_num_variables, &folding_factor),
        folding_factor,
        soundness_type: SecurityAssumption::CapacityBound,
        starting_log_inv_rate: 1,
    };
    let whir_config = WhirConfig::new(stacked_num_variables, params).unwrap();
    MultiFibConfig {
        pcs: MyPcs::new(whir_config, MyDft::default(), mmcs),
    }
}

/// A challenger seeded with the same domain separator on both prove and verify sides.
pub fn challenger(config: &MultiFibConfig) -> MyChallenger {
    let mut challenger = MyChallenger::new(perm());
    let mut ds = DomainSeparator::new(vec![]);
    config.pcs.add_domain_separator::<8>(&mut ds);
    ds.observe_domain_separator(&mut challenger);
    challenger
}

// ---------------------------------------------------------------------------
// M1: the same multilinear AIR path, but committed with the trivial
// field-agnostic PCS instead of WHIR. Nothing above WHIR (the AIR, the driver)
// changes — only the PCS type and its (WHIR-free) challenger setup.
// ---------------------------------------------------------------------------

/// The trivial PCS reuses the same BabyBear Merkle scheme purely to produce an
/// observable commitment; its opening path is field-agnostic (no two-adicity).
pub type TrivialPcs = TrivialMultilinearPcs<F, EF, MyMmcs>;

/// A trivial-PCS-backed multilinear AIR configuration over BabyBear.
pub struct MultiFibConfigTrivial {
    pcs: TrivialPcs,
}

impl MultiStarkConfig for MultiFibConfigTrivial {
    type Val = F;
    type Challenge = EF;
    type Challenger = MyChallenger;
    type Pcs = TrivialPcs;

    fn pcs(&self) -> &TrivialPcs {
        &self.pcs
    }

    fn min_num_variables(&self) -> usize {
        // The trivial PCS commits tables at their natural arity (no padding),
        // so any positive arity is acceptable.
        1
    }

    fn build_witness(&self, tables: Vec<Table<F>>) -> Vec<Table<F>> {
        // The trivial PCS's witness is just the committed tables, unstacked.
        tables
    }

    fn committed_table<'a>(
        &self,
        prover_data: &'a Vec<Table<F>>,
        table_index: usize,
    ) -> &'a Table<F> {
        &prover_data[table_index]
    }
}

/// Build a trivial-PCS configuration sized for a single trace.
pub fn config_for_trivial(log_height: usize, width: usize) -> MultiFibConfigTrivial {
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let num_vars = log_height + log2_ceil_usize(width);
    MultiFibConfigTrivial {
        pcs: TrivialPcs::new(mmcs, num_vars),
    }
}

/// The trivial PCS injects no domain separator, so a bare seeded challenger
/// suffices — identical on prove and verify sides.
pub fn challenger_trivial() -> MyChallenger {
    MyChallenger::new(perm())
}

// ---------------------------------------------------------------------------
// M4: the same multilinear AIR path, committed with the **real recursive
// Brakedown** tensor PCS (linear-time encoding — genuinely Brakedown, not the
// RS/Ligero base). The stacked-layout sumcheck reduction is shared with WHIR;
// only the final `P(r) = v` step is discharged by the tensor opening (no FFT,
// but the reduction here happens to reuse two-adicity via `PrefixProver`).
// ---------------------------------------------------------------------------

/// The real recursive Brakedown linear-time code the tensor PCS commits through.
pub type BrakedownCode = RealBrakedownCode<F, EF>;
/// The standalone tensor commitment wrapped by the multilinear adapter.
pub type BrakedownTensor = TensorPcs<F, EF, BrakedownCode, MyMmcs>;
/// The Brakedown-backed multilinear PCS used by multi-stark.
pub type BrakedownPcs = BrakedownMlPcs<F, EF, BrakedownCode, MyMmcs, L>;

/// Fixed seed for the Brakedown code's random `A`/`B` expander matrices.
///
/// `config_for_brakedown` is called independently on the prover and verifier
/// sides; seeding the code RNG deterministically guarantees both build the
/// **identical** linear code (same sparse blocks), so their commitments,
/// encodings, and query checks line up.
const BRAKEDOWN_CODE_SEED: u64 = 0xB2A4E_D0;

/// Preprocessing/padding floor for the Brakedown config. The tensor PCS commits
/// the stacked polynomial at its natural arity, so any positive arity is fine.
const BRAKEDOWN_FOLDING: usize = 1;

/// A Brakedown-PCS-backed multilinear AIR configuration over BabyBear.
pub struct MultiFibConfigBrakedown {
    pcs: BrakedownPcs,
}

impl MultiStarkConfig for MultiFibConfigBrakedown {
    type Val = F;
    type Challenge = EF;
    type Challenger = MyChallenger;
    type Pcs = BrakedownPcs;

    fn pcs(&self) -> &BrakedownPcs {
        &self.pcs
    }

    fn min_num_variables(&self) -> usize {
        BRAKEDOWN_FOLDING
    }

    fn build_witness(&self, tables: Vec<Table<F>>) -> Witness<F> {
        L::new_witness(tables, BRAKEDOWN_FOLDING)
    }

    fn committed_table<'a>(
        &self,
        prover_data: &'a p3_brakedown::BrakedownProverData<F, EF, MyMmcs, L>,
        table_index: usize,
    ) -> &'a Table<F> {
        prover_data.table(table_index)
    }
}

/// Build a Brakedown-PCS configuration sized for a single trace.
///
/// The committed polynomial stacks every trace column, so the tensor scheme is
/// sized for the stacked variable count `log_height + ceil(log2(width))`, split
/// into `log_rows = stacked / 2` row variables and the rest as columns.
pub fn config_for_brakedown(log_height: usize, width: usize) -> MultiFibConfigBrakedown {
    let stacked_num_variables = log_height + log2_ceil_usize(width);
    let log_rows = stacked_num_variables / 2;
    // Columns carry the remaining variables; a power of two so the row/column
    // split is itself a clean multilinear tensor.
    let cols = 1 << (stacked_num_variables - log_rows);

    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    // Real recursive Brakedown code over the (power-of-two) message length `cols`.
    // Built from a fixed seed so the prover and verifier configs share one code.
    let mut code_rng = SmallRng::seed_from_u64(BRAKEDOWN_CODE_SEED);
    let code = BrakedownCode::from_message_len(cols, &mut code_rng);
    // Brakedown's codeword length is not a power of two; clamp the query count to it.
    let num_queries = core::cmp::min(32, code.codeword_len());
    let tensor = TensorPcs::new(code, mmcs, log_rows, num_queries);
    MultiFibConfigBrakedown {
        pcs: BrakedownPcs::new(tensor),
    }
}

/// The Brakedown PCS injects no domain separator, so a bare seeded challenger
/// suffices — identical on prove and verify sides.
pub fn challenger_brakedown() -> MyChallenger {
    MyChallenger::new(perm())
}
