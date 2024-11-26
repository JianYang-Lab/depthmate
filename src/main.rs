use anyhow::Result;
use cli::parse_cli;
use csv::Writer;
use depth::{DefaultReadFilter, DepthProcessor, PilePosition, ReadFilter};
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

    let read_filter = DefaultReadFilter::new(args.include_flags, args.exclude_flags, args.min_mapq);

    let tumor_bam_path = PathBuf::from(args.tumor);
    let normal_bam_path = PathBuf::from(args.normal);

    let depth_processer = DepthProcessor::new(tumor_bam_path, normal_bam_path, read_filter);

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    let final_pile_res = variants
        .par_iter()
        .flat_map(|variant| {
            get_pile_res(variant, &depth_processer).unwrap_or_else(|e| {
                eprintln!("Error: {:?}", e);
                vec![]
            })
        })
        .collect::<Vec<PilePosition>>();

    let mut writer = Writer::from_path(args.output)?;
    for pile in final_pile_res {
        writer.serialize(pile)?;
    }

    Ok(())
}

fn get_pile_res<F>(
    variant: &Variant,
    depth_processer: &DepthProcessor<F>,
) -> Result<Vec<PilePosition>>
where
    F: ReadFilter + Send,
{
    let chr = variant.chrom.clone();
    let pos = variant.pos as u32;
    let ref_base = variant.ref_base.as_bytes()[0];
    let alt_base = variant.alt_base.as_bytes()[0];

    let res = depth_processer.process(&chr, pos, ref_base, alt_base)?;
    Ok(res)
}
