/// Sequence based manipulations on a GFF3/fasta combination.
pub mod seq;
/// Statistical manipulations on a GFF3/fasta combination.
pub mod stat;
/// Utility functions used in other modules.
pub mod utils;

/// For easier access to degeneracy enum
pub use utils::Degeneracy;