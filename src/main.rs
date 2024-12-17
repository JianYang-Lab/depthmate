use anyhow::Result;
use cli::parse_cli;
use csv::Writer;
use depth::{DefaultReadFilter, DepthProcessor, PilePosition, ReadFilter};
use itertools::Itertools;
use rayon::prelude::*;
use std::path::PathBuf;
use vcf::{get_variant, Variant};
mod cli;
mod depth;
mod vcf;

fn main() -> Result<()> {
    let args = parse_cli()?;

    let vcf_path = args.vcf;
    let variants = get_variant(&vcf_path)?;
    println!("Number of variants: {}", variants.len());

    let read_filter = DefaultReadFilter::new(args.include_flags, args.exclude_flags, args.min_mapq);

    let tumor_bam_path = PathBuf::from(args.tumor);
    let normal_bam_path = PathBuf::from(args.normal);

    let depth_processer = DepthProcessor::new(tumor_bam_path, normal_bam_path, read_filter);

    let t_n_types = [78, 84];
    // product of variants and t_n_types
    // vars: [v1,v2,v3] t_n_types: [78, 84]
    // Cartesian product and parallel process
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    let final_pile_res = variants
        .iter()
        .cartesian_product(t_n_types.iter())
        .par_bridge()
        .map(|(variant, tntype)| get_pile_res(variant, &depth_processer, *tntype))
        .collect::<Result<Vec<Vec<PilePosition>>>>()?;

    // get result from two thread pool

    let mut writer = Writer::from_path(args.output)?;
    for pile_vec in final_pile_res {
        for pile in pile_vec {
            writer.serialize(pile)?;
        }
    }

    Ok(())
}

fn get_pile_res<F>(
    variant: &Variant,
    depth_processer: &DepthProcessor<F>,
    tntype: u8,
) -> Result<Vec<PilePosition>>
where
    F: ReadFilter + Send,
{
    let chr = variant.chrom.clone();
    let pos = variant.pos as u32;
    let ref_base = variant.ref_base.as_bytes()[0];
    let alt_base = variant.alt_base.as_bytes()[0];

    let res = depth_processer.process(&chr, pos, ref_base, alt_base, tntype)?;
    Ok(res)
}
