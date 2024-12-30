//! This feature is inspired by [perbase](https://github.com/sstadick/perbase)
use anyhow::{Ok, Result};
use rust_htslib::bam::pileup::{Alignment, Pileup};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::HeaderView;
use rust_htslib::{bam, bam::Read};
use serde::{Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use std::default;
use std::path::PathBuf;

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct PilePosition {
    pub chrom: String,
    pub pos: usize,
    pub all_depth: usize,
    pub alt_depth: usize,
    pub ref_depth: usize,
    pub other_depth: usize,
    pub indel_n_lowq: usize,
    pub tntype: String,
    pub vaf: f64,
    pub barcode: String,
    // pub umi: String,
}

pub trait Position: Default {
    /// Create a new position with all other values zeroed
    fn new(ref_seq: String, pos: u32) -> Self;
}

impl Position for PilePosition {
    fn new(ref_seq: String, pos: u32) -> Self {
        Self {
            chrom: ref_seq,
            pos: pos as usize,
            ..default::Default::default()
        }
    }
}

impl PilePosition {
    /// Given a record, update the counts at this position
    #[inline(always)]
    fn update<F: ReadFilter>(
        &mut self,
        alignment: &Alignment,
        record: Record,
        read_filter: &F,
        ref_base: u8,
        alt_base: u8,
    ) {
        // // get UMI and barcode
        // // barcode tag: CB
        // // UMI tag: UB

        let barcode = match record.aux(b"CB") {
            Result::Ok(bam::record::Aux::String(b)) => b.to_string(),
            _ => "NA".to_string(),
        };

        // let umi = match record.aux(b"UB") {
        //     ResultOk(bam::record::Aux::String(b)) => b.to_string(),
        //     _ => "NA".to_string(),
        // };

        // self.barcode.push_str(&barcode);
        // self.umi.push_str(&umi);

        // self.barcode.push(',');
        // self.umi.push(',');

        if !read_filter.filter_read(&record, Some(alignment)) {
            self.all_depth -= 1;
            self.indel_n_lowq += 1;
            return;
        }
        // NB: Order matters here, a refskip is true for both is_del and is_refskip
        // while a true del is only true for is_del
        if alignment.is_refskip() | alignment.is_del() {
            self.indel_n_lowq += 1;
            self.all_depth -= 1;
        } else {
            // We have an actual base!

            // check ref_base and lat_base

            let curr_base = record.seq()[alignment.qpos().unwrap()];

            if curr_base == ref_base {
                self.ref_depth += 1;
            } else if curr_base == alt_base {
                self.alt_depth += 1;
                self.barcode.push_str(&barcode);
                self.barcode.push(',');
            } else {
                self.other_depth += 1;
            }

            // Check for insertions
            if let bam::pileup::Indel::Ins(_len) = alignment.indel() {
                self.indel_n_lowq += 1;
            }
        }

        // calculate VAF
        self.vaf = self.alt_depth as f64 / self.all_depth as f64;
    }

    #[inline]
    pub fn from_pileup<F: ReadFilter>(
        pileup: Pileup,
        header: &bam::HeaderView,
        read_filter: &F,
        ref_base: u8,
        alt_base: u8,
        tntype: u8,
    ) -> Self {
        let name = Self::compact_refseq(header, pileup.tid());
        // make output 1-based
        let mut pos = Self::new(name.to_string(), pileup.pos() + 1);
        pos.all_depth = pileup.depth() as usize;
        pos.tntype = match tntype {
            84 => "T".to_string(),
            78 => "N".to_string(),
            _ => "U".to_string(),
        };

        for alignment in pileup.alignments() {
            let record = alignment.record();
            Self::update(
                &mut pos,
                &alignment,
                record,
                read_filter,
                ref_base,
                alt_base,
            );
        }
        pos
    }

    #[inline]
    pub fn compact_refseq(header: &HeaderView, tid: u32) -> SmartString<LazyCompact> {
        let name = std::str::from_utf8(header.tid2name(tid)).unwrap();
        smartstring::alias::String::from(name)
    }
}

/// Anything that implements ReadFilter can apply a filter set to read.
pub trait ReadFilter {
    /// filters a read, true is pass, false if fail
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool;
}

/// A straightforward read filter.
pub struct DefaultReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

impl DefaultReadFilter {
    /// Create an OnlyDepthReadFilter
    pub fn new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self {
        Self {
            include_flags,
            exclude_flags,
            min_mapq,
        }
    }
}

impl ReadFilter for DefaultReadFilter {
    /// Filter reads based SAM flags and mapping quality
    #[inline(always)]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        let flags = read.flags();
        (!flags) & self.include_flags == 0
            && flags & self.exclude_flags == 0
            && read.mapq() >= self.min_mapq
    }
}

pub(crate) struct DepthProcessor<F: ReadFilter + Send> {
    /// path to indexed tumor BAM
    pub tumor_reads: PathBuf,
    /// path to indexed normal BAM
    pub normal_reads: PathBuf,
    /// implementation of [position::ReadFilter] that will be used
    pub read_filter: F,
}

impl<F: ReadFilter + Send> DepthProcessor<F> {
    /// Create a new OnlyDepthProcessor
    pub fn new(tumor_reads: PathBuf, normal_reads: PathBuf, read_filter: F) -> Self {
        Self {
            tumor_reads,
            normal_reads,
            read_filter,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn process_region(
        &self,
        chrom: &str,
        pos: u32,
        bam_path: PathBuf,
        ref_base: u8,
        alt_base: u8,
        tntype: u8, // T or N
        res: &mut Vec<PilePosition>,
    ) -> Result<()> {
        // Create a reader
        let mut reader = bam::IndexedReader::from_path(bam_path)?;
        let header = reader.header().to_owned();
        // fetch the region of interest
        reader.fetch((chrom, pos - 1, pos))?;

        let pileup = reader.pileup();

        pileup.for_each(|p| {
            let pileup = p.expect("msg");
            if pileup.pos() + 1 == pos {
                let pos = PilePosition::from_pileup(
                    pileup,
                    &header,
                    &self.read_filter,
                    ref_base,
                    alt_base,
                    tntype,
                );
                res.push(pos);
            }
        });
        Ok(())
    }

    /// process tumor and normal reads for a given region
    pub fn process(
        &self,
        chrom: &str,
        pos: u32,
        ref_base: u8,
        alt_base: u8,
        tntype: u8,
    ) -> Result<Vec<PilePosition>> {
        let mut output = Vec::new();
        match tntype {
            84 => self.process_region(
                chrom,
                pos,
                self.tumor_reads.clone(),
                ref_base,
                alt_base,
                tntype,
                &mut output,
            )?,
            78 => self.process_region(
                chrom,
                pos,
                self.normal_reads.clone(),
                ref_base,
                alt_base,
                tntype,
                &mut output,
            )?,
            _ => {
                anyhow::bail!("Invalid tntype");
            }
        }

        // self.process_region(
        //     chrom,
        //     pos,
        //     self.normal_reads.clone(),
        //     ref_base,
        //     alt_base,
        //     78,
        //     &mut output,
        // )?;

        Ok(output)
    }
}
