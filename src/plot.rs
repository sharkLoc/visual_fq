use plotters::prelude::*;
use std::collections::HashMap;

#[allow(non_snake_case)]
pub fn plot_base(
    data: Vec<HashMap<u64, f64>>,
    max_readlen: u64,
    figBase: String,
    width: u32,
    height: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    // line plot for base A T G C N rate in position
    let n = format!("{}.png", figBase);
    let root = BitMapBackend::new(n.as_str(), (width, height)).into_drawing_area();

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

    let (mut nt_a, mut nt_t, mut nt_g, mut nt_c, mut nt_n) =
        (vec![], vec![], vec![], vec![], vec![]);
    for k in 1..=max_readlen {
        nt_a.push((k as f32, data[0][&k] as f32));
        nt_t.push((k as f32, data[1][&k] as f32));
        nt_g.push((k as f32, data[2][&k] as f32));
        nt_c.push((k as f32, data[3][&k] as f32));
        nt_n.push((k as f32, data[4][&k] as f32));
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

    Ok(())
}

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
pub fn plot_qual(
    data: Vec<Vec<u64>>,
    max_qvalue: u8,
    max_readlen: u64,
    figQual: String,
    width: u32,
    height: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    // heatmap plot for base qualtity
    let mut plot_heatmap: Vec<(u64, u8, u64)> = Vec::new();
    let mut max = 0;
    let mut row = 0usize;
    for x in data {
        row += 1;
        for j in 0..=max_qvalue {
            let cot = x[j as usize];
            if max < cot {
                max = cot;
            }
            plot_heatmap.push((row as u64, j, cot));
        }
    }

    let n = format!("{}.png", figQual);
    let heatmap = BitMapBackend::new(n.as_str(), (width, height)).into_drawing_area();
    heatmap.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&heatmap)
        .caption("quality distrbution plot", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0u64..(max_readlen + 1), 0u64..(max_qvalue + 1) as u64)?;

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
            .map(move |(x, y, z)| (*x, *y as u64, *z as f64))
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
