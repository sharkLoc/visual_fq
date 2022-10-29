use clap::Parser;

///title
#[allow(non_snake_case)]
#[derive(Parser, Debug)]
#[clap(setting = clap::AppSettings::DeriveDisplayOrder)]
pub struct Opt {
    ///input fastq1 file
    #[clap(short = 'i', long)]
    pub read1: String,

    ///input fastq2 file [for PE reads]
    #[clap(short = 'I', long)]
    pub read2: Option<String>,

    ///output detail matrix file
    #[clap(short = 'o', long)]
    pub out: String,

    ///output summary log file
    #[clap(short = 'O', long)]
    pub log: String,

    ///output base figure prefix name
    #[clap(short = 'b', long, default_value = "base_plot")]
    pub figBase: String,

    ///output quality figure prefix name
    #[clap(short = 'q', long, default_value = "qual_plot")]
    pub figQual: String,

    ///set output figure width
    #[clap(short = 'w', long, default_value = "960")]
    pub width: u32,

    ///set output figure height
    #[clap(short = 'h', long, default_value = "540")]
    pub height: u32,
}
