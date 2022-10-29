mod cli;
mod data;
mod plot;

use clap::Parser;
use cli::Opt;
use data::read_fq;

#[allow(unused_variables)]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opt = cli::Opt::parse();
    match opt {
        Opt {
            read1,
            read2,
            out,
            log,
            figBase,
            figQual,
            width,
            height,
        } => {
            read_fq(read1, read2, out, figBase, figQual, width, height)?;
        }
    }

    Ok(())
}
