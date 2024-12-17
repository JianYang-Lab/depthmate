use anyhow::Result;
use clap::{command, Parser};

#[derive(Parser, Debug)]
#[command(name = "depthmate")]
#[command(about = "tumor-normal base depth")]
#[command(long_about = "long_about todo!!!")]
#[command(author, version)]
#[command(
    help_template = "{name} -- {about}\n\nVersion: {version}\n\nAuthors: {author}\
    \n\n{usage-heading} {usage}\n\n{all-args}"
)]
pub struct Cli {
    /// Input normal bam file with index
    #[arg(short, long, help_heading = Some("I/O Options"))]
    pub normal: String,
    /// Input tumor bam file with index
    #[arg(short, long, help_heading = Some("I/O Options"))]
    pub tumor: String,
    /// input vcf file
    #[arg(short, long, help_heading = Some("I/O Options"))]
    pub vcf: String,
    /// Output file
    #[arg(short, long, help_heading = Some("I/O Options"))]
    pub output: String,

    /// Included flags
    #[arg(default_value = "0", short, long, help_heading = Some("Filter Options"))]
    pub include_flags: u16,
    /// Excluded flags
    #[arg(default_value = "0", short, long, help_heading = Some("Filter Options"))]
    pub exclude_flags: u16,
    /// Minimum mapping quality
    #[arg(default_value = "20", short, long, help_heading = Some("Filter Options"))]
    pub min_mapq: u8,

    /// Threads
    #[arg(default_value = "1", short = '@', long)]
    pub threads: usize,
}

pub fn parse_cli() -> Result<Cli> {
    let cli = Cli::parse();
    Ok(cli)
}
