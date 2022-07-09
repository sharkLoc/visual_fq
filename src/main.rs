use plotters::prelude::*;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() == 1 {
        println!("[Usage information]\nuncompressed fastq:\n\t{} your_file.fq > out.matrix 2> summary.log\n",&args[0]);
        println!("compressed fastq:\n\t gzip -dc your_file.fq.gz | {} /dev/stdin > out.matrix 2> summary.log\n",&args[0]);
        println!("result files:\n\t out.matrix  summary.log  Base_plot.png Qual_plot.png\n");
        std::process::exit(1);
    }

    let fp = File::open(&args[1])?;
    let reader = BufReader::new(&fp);

    let mut max_readlen = 0u32;
    let mut max_qvalue = 0u8;
    let mut pos_info: HashMap<u32, HashMap<u8, u32>> = HashMap::new();
    let mut pos: HashMap<u8, u32> = HashMap::new();
    let mut nt_info: HashMap<u32, HashMap<char, u32>> = HashMap::new();
    let mut nt: HashMap<char, u32> = HashMap::new();
    let mut total_seq = 0u32;
    let mut total_nt = 0u32;
    let mut total_gc = 0u32;

    let mut n = 0u8;
    for line in reader.lines() {
        n += 1;
        if n == 1 || n == 3 {
            continue;
        }
        if n == 2 {
            if let Ok(seq) = &line {
                total_seq += 1;
                for (order, base) in seq.as_bytes().iter().enumerate() {
                    let real_p = order as u32 + 1;
                    let base = *base as char;
                    *nt.entry(base).or_insert(0) += 1;

                    if nt_info.contains_key(&real_p) {
                        *nt_info.get_mut(&real_p).unwrap().entry(base).or_insert(0) += 1;
                    } else {
                        nt_info.insert(real_p, nt.clone());
                    }
                    nt.clear();
                }
            }
        }
        if n == 4 {
            if let Ok(qual) = &line {
                let length = qual.len() as u32;
                if length > max_readlen {
                    max_readlen = length;
                }

                for (order, qnum) in qual.as_bytes().iter().enumerate() {
                    let q = qnum - 33;
                    if q > max_qvalue {
                        max_qvalue = q;
                    }

                    let real_p = order as u32 + 1;
                    *pos.entry(q).or_insert(0) += 1;

                    if pos_info.contains_key(&real_p) {
                        *pos_info.get_mut(&real_p).unwrap().entry(q).or_insert(0) += 1;
                    } else {
                        pos_info.entry(real_p).or_insert(pos.clone());
                    }
                    pos.clear();
                }
            }
            n = 0;
        }
    }

    //for head
    for m in 0..=max_qvalue {
        if m == 0 {
            print!("Iterm\tTotal\tA\tT\tG\tC\tN\t");
        }
        if m == max_qvalue {
            print!("{}\n", &m);
        } else {
            print!("{}\t", &m);
        }
    }
    let mut plot_base_a: HashMap<u32, f64> = HashMap::new();
    let mut plot_base_t: HashMap<u32, f64> = HashMap::new();
    let mut plot_base_g: HashMap<u32, f64> = HashMap::new();
    let mut plot_base_c: HashMap<u32, f64> = HashMap::new();
    let mut plot_base_n: HashMap<u32, f64> = HashMap::new();

    for i in 0..max_readlen {
        let x = i + 1;
        if let Some(tmp) = pos_info.get(&x) {
            if let Some(tmp2) = nt_info.get(&x) {
                let base_a = match tmp2.get(&'A') {
                    Some(n) => n,
                    None => &0,
                };
                let base_t = match tmp2.get(&'T') {
                    Some(n) => n,
                    None => &0,
                };
                let base_g = match tmp2.get(&'G') {
                    Some(n) => n,
                    None => &0,
                };
                let base_c = match tmp2.get(&'C') {
                    Some(n) => n,
                    None => &0,
                };
                let base_n = match tmp2.get(&'N') {
                    Some(n) => n,
                    None => &0,
                };
                let total = base_a + base_t + base_g + base_c + base_n;
                total_nt += total;
                total_gc += base_g + base_c;
                let rate_a = *base_a as f64 / total as f64 * 100.0;
                let rate_t = *base_t as f64 / total as f64 * 100.0;
                let rate_g = *base_g as f64 / total as f64 * 100.0;
                let rate_c = *base_c as f64 / total as f64 * 100.0;
                let rate_n = *base_n as f64 / total as f64 * 100.0;
                plot_base_a.insert(x, rate_a);
                plot_base_t.insert(x, rate_t);
                plot_base_g.insert(x, rate_g);
                plot_base_c.insert(x, rate_c);
                plot_base_n.insert(x, rate_n);
                print!(
                    "pos:{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t",
                    &x, total, rate_a, rate_t, rate_g, rate_c, rate_n
                );
            }
            for j in 0..=max_qvalue {
                match tmp.get(&j) {
                    Some(x) => {
                        if j == max_qvalue {
                            print!("{}\n", x);
                        } else {
                            print!("{}\t", x);
                        }
                    }
                    None => {
                        if j == max_qvalue {
                            print!("{}\n", 0);
                        } else {
                            print!("{}\t", 0);
                        }
                    }
                }
            }
        }
    }

    eprintln!("Total read number:\t{}", total_seq);
    eprintln!("Total base number:\t{}", total_nt);
    eprintln!("Max read length:\t{}", max_readlen);
    eprintln!(
        "Reads average length:\t{:.1}",
        total_nt as f64 / total_seq as f64
    );
    eprintln!(
        "Total GC content:\t{:.2}",
        total_gc as f64 / total_nt as f64 * 100.0
    );

    // line plot for base A T G C N rate in position
    let root = BitMapBackend::new("Base_plot.png", (960, 540)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);

    let mut chart = ChartBuilder::on(&root)
        .caption("Base distrbution plot", ("sans-serif", 40).into_font())
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-0.1f32..(max_readlen + 1) as f32, -0.5f32..51f32)?;

    chart
        .configure_mesh()
        .x_labels(20)
        .x_desc("position")
        .x_label_formatter(&|x| format!("{:.0}", x))
        .y_labels(10)
        .y_label_formatter(&|x| format!("{:.1}", x))
        .y_desc("percent")
        .draw()?;

    let mut nt_a = vec![];
    let mut nt_t = vec![];
    let mut nt_g = vec![];
    let mut nt_c = vec![];
    let mut nt_n = vec![];
    for k in 1..=max_readlen {
        nt_a.push((k as f32, plot_base_a[&k] as f32));
        nt_t.push((k as f32, plot_base_t[&k] as f32));
        nt_g.push((k as f32, plot_base_g[&k] as f32));
        nt_c.push((k as f32, plot_base_c[&k] as f32));
        nt_n.push((k as f32, plot_base_n[&k] as f32));
    }

    chart
        .draw_series(LineSeries::new(nt_a, &RED))
        .unwrap()
        .label("A")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    chart
        .draw_series(LineSeries::new(nt_t, &GREEN))
        .unwrap()
        .label("T")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
    chart
        .draw_series(LineSeries::new(nt_g, &YELLOW))
        .unwrap()
        .label("G")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &YELLOW));
    chart
        .draw_series(LineSeries::new(nt_c, &BLACK))
        .unwrap()
        .label("C")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));
    chart
        .draw_series(LineSeries::new(nt_n, &BLUE))
        .unwrap()
        .label("N")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.9))
        .border_style(&BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    // heatmap plot for base qualtity
    let mut plot_heatmap: Vec<(u32, u8, u32)> = Vec::new();
    let mut max = 0;
    for i in 0..max_readlen {
        let x = i + 1;
        if let Some(tmp) = pos_info.get(&x) {
            for j in 0..=max_qvalue {
                let count = match tmp.get(&j) {
                    Some(v) => v,
                    None => &0,
                };
                if max < *count {
                    max = *count;
                }
                plot_heatmap.push((x, j, *count));
            }
        }
    }
    let heatmap = BitMapBackend::new("Qual_plot.png", (960, 540)).into_drawing_area();
    heatmap.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&heatmap)
        .caption("quality distrbution plot", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0u32..(max_readlen + 1), 0u32..(max_qvalue + 1) as u32)?;

    chart
        .configure_mesh()
        .x_desc("position")
        .x_labels(20)
        .y_labels(20)
        .x_label_offset(0)
        .y_label_offset(0)
        .disable_x_mesh()
        .disable_y_mesh()
        .label_style(("sans-serif", 20))
        .draw()?;

    chart.draw_series(
        plot_heatmap
            .iter()
            .map(move |(x, y, z)| (*x, *y as u32, *z as f64))
            .map(|(x, y, z)| {
                Rectangle::new(
                    [(x, y), (x + 1, y + 1)],
                    HSLColor(
                        20.0 / 360.0 * (1.0 - z / max as f64), // 单色系渐变
                        1.0,                                   // 灰度（对比度）
                        1.0 - z / max as f64, // 亮度 (明暗度: black :0 --> white: 1)
                    )
                    .filled(),
                )
            }),
    )?;

    Ok(())
}
