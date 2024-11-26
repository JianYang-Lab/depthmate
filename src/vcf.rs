use anyhow::{Ok, Result};
use noodles::vcf::{self, variant::record::AlternateBases};

// get chrom,pos,ref_base,alt_base from a vcf path

#[derive(Debug)]
pub struct Variant {
    pub chrom: String,
    pub pos: usize,
    pub ref_base: String,
    pub alt_base: String,
}

pub fn get_variant(vcf_path: &str) -> Result<Vec<Variant>> {
    let mut variants = Vec::new();
    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path)?;
    let _header = reader.read_header()?;

    for rec in reader.records() {
        let rec = rec?;
        let pos = rec.variant_start().transpose()?.unwrap().get();
        let chrom = rec.reference_sequence_name().to_string();
        let ref_base = rec.reference_bases().to_string();
        let alt_base = rec.alternate_bases();
        let alt_base = alt_base.iter().next().unwrap()?.to_string();

        let variant = Variant {
            chrom,
            pos,
            ref_base,
            alt_base,
        };
        variants.push(variant);
    }

    Ok(variants)
}
