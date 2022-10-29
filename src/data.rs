use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use crate::plot::{plot_base, plot_qual};

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
type hash_base = HashMap<char, u64>;
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
type hash_base_pos = HashMap<u64, hash_base>;
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
type hash_qual = HashMap<u8, u64>;
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
type hash_qual_pos = HashMap<u64, hash_qual>;

pub type Plotline = HashMap<u64, f64>;

#[allow(non_snake_case)]
pub fn read_fq(
    r1: String,
    r2: Option<String>,
    out: String,
    figBase: String,
    figQual: String,
    width: u32,
    height: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    match r2 {
        Some(r2) => {
            open_pe(r1, r2, out, figBase, figQual, width, height)?;
        }
        None => {
            open_se(r1, out, figBase, figQual, width, height)?;
        }
    }
    Ok(())
}

#[allow(unused_variables)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
pub fn open_pe(
    r1: String,
    r2: String,
    out: String,
    figBase: String,
    figQual: String,
    width: u32,
    height: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let fp1 = File::open(&r1)?;
    let fp2 = File::open(&r2)?;
    let reader1 = BufReader::new(&fp1);
    let reader2 = BufReader::new(&fp2);

    let mut nt_info1: hash_base_pos = HashMap::new();
    let mut nt_info2: hash_base_pos = HashMap::new();
    let mut nt1: hash_base = HashMap::new();
    let mut nt2: hash_base = HashMap::new();

    let mut pos_info1: hash_qual_pos = HashMap::new();
    let mut pos_info2: hash_qual_pos = HashMap::new();
    let mut pos1: hash_qual = HashMap::new();
    let mut pos2: hash_qual = HashMap::new();

    let (mut n, mut max_qvalue1, mut max_qvalue2) = (0u8, 0u8, 0u8);
    let (mut max_readlen1, mut max_readlen2) = (0u64, 0u64);
    let (mut total_seq1, total_nt1, total_gc1) = (0u64, 0u64, 0u64);
    let (mut total_seq2, total_nt2, total_gc2) = (0u64, 0u64, 0u64);
    for ret in reader1.lines().zip(reader2.lines()) {
        n += 1;
        if n == 1 || n == 3 {
            continue;
        }
        if n == 2 {
            if let Ok(s1) = &ret.0 {
                total_seq1 += 1;
                for (order, base) in s1.as_bytes().iter().enumerate() {
                    let real_p = order as u64 + 1;
                    let base = *base as char;
                    *nt1.entry(base).or_insert(0) += 1;

                    if nt_info1.contains_key(&real_p) {
                        *nt_info1.get_mut(&real_p).unwrap().entry(base).or_insert(0) += 1;
                    } else {
                        nt_info1.insert(real_p, nt1.clone());
                    }
                    nt1.clear();
                }
                if let Ok(s2) = &ret.1 {
                    total_seq2 += 1;
                    for (order, base) in s2.as_bytes().iter().enumerate() {
                        let real_p = order as u64 + 1;
                        let base = *base as char;
                        *nt2.entry(base).or_insert(0) += 1;

                        if nt_info2.contains_key(&real_p) {
                            *nt_info2.get_mut(&real_p).unwrap().entry(base).or_insert(0) += 1;
                        } else {
                            nt_info2.insert(real_p, nt2.clone());
                        }
                        nt2.clear();
                    }
                }
            }
        }
        if n == 4 {
            if let Ok(qual) = &ret.0 {
                let length = qual.len() as u64;
                if length > max_readlen1 {
                    max_readlen1 = length;
                }

                for (order, qnum) in qual.as_bytes().iter().enumerate() {
                    let q = qnum - 33;
                    if q > max_qvalue1 {
                        max_qvalue1 = q;
                    }

                    let real_p = order as u64 + 1;
                    *pos1.entry(q).or_insert(0) += 1;

                    if pos_info1.contains_key(&real_p) {
                        *pos_info1.get_mut(&real_p).unwrap().entry(q).or_insert(0) += 1;
                    } else {
                        pos_info1.entry(real_p).or_insert(pos1.clone());
                    }
                    pos1.clear();
                }
            }
            if let Ok(qual) = &ret.1 {
                let length = qual.len() as u64;
                if length > max_readlen2 {
                    max_readlen2 = length;
                }

                for (order, qnum) in qual.as_bytes().iter().enumerate() {
                    let q = qnum - 33;
                    if q > max_qvalue2 {
                        max_qvalue2 = q;
                    }

                    let real_p = order as u64 + 1;
                    *pos2.entry(q).or_insert(0) += 1;

                    if pos_info2.contains_key(&real_p) {
                        *pos_info2.get_mut(&real_p).unwrap().entry(q).or_insert(0) += 1;
                    } else {
                        pos_info2.entry(real_p).or_insert(pos2.clone());
                    }
                    pos2.clear();
                }
            }
            n = 0;
        }
    }

    // plot data
    let max_qvalue: u8;
    let mut qual = vec![];
    if max_qvalue1 > max_qvalue2 {
        max_qvalue = max_qvalue1;
    } else {
        max_qvalue = max_qvalue2;
    }
    let qual1 = format_qual(&pos_info1, max_readlen1, max_qvalue1)?;
    let qual2 = format_qual(&pos_info2, max_readlen2, max_qvalue2)?;
    for x in qual1 {
        qual.push(x);
    }
    for x in qual2 {
        qual.push(x);
    }
    plot_qual(
        qual.clone(),
        max_qvalue,
        max_readlen1 + max_readlen2,
        figQual,
        width * 2,
        height,
    )?;

    let base1 = format_base(&nt_info1, max_readlen1)?;
    let base2 = format_base(&nt_info2, max_readlen2)?;
    let mut base: Vec<HashMap<u64, f64>> = vec![];
    for i in 0..5 {
        base.push(base1[i].clone());
        for (k, v) in base2[i].clone().into_iter() {
            base[i].insert(k + max_readlen1, v);
        }
    }
    plot_base(
        base.clone(),
        max_readlen1 + max_readlen2,
        figBase,
        width * 2,
        height,
    )?;
    show_mat(out, max_qvalue, max_readlen1 + max_readlen2, base, qual)?;

    Ok(())
}

fn format_qual(
    x: &hash_qual_pos,
    max_readlen: u64,
    max_qvalue: u8,
) -> Result<Vec<Vec<u64>>, Box<dyn std::error::Error>> {
    let mut vec_qual_pos = vec![];
    for i in 0..max_readlen {
        let pos = i + 1;
        let mut vec_qual = vec![];
        if let Some(qu) = x.get(&pos) {
            for j in 0..=max_qvalue {
                let cot = match qu.get(&j) {
                    Some(x) => x,
                    None => &0,
                };
                vec_qual.push(*cot);
            }
        }
        vec_qual_pos.push(vec_qual);
    }

    Ok(vec_qual_pos)
}

fn format_base(
    y: &hash_base_pos,
    max_readlen: u64,
) -> Result<Vec<HashMap<u64, f64>>, Box<dyn std::error::Error>> {
    let mut plot_base_a: Plotline = HashMap::new();
    let mut plot_base_t: Plotline = HashMap::new();
    let mut plot_base_g: Plotline = HashMap::new();
    let mut plot_base_c: Plotline = HashMap::new();
    let mut plot_base_n: Plotline = HashMap::new();
    for i in 0..max_readlen {
        let pos = i + 1;
        if let Some(cot) = y.get(&pos) {
            let base_a = match cot.get(&'A') {
                Some(n) => n,
                None => &0,
            };
            let base_t = match cot.get(&'T') {
                Some(n) => n,
                None => &0,
            };
            let base_g = match cot.get(&'G') {
                Some(n) => n,
                None => &0,
            };
            let base_c = match cot.get(&'C') {
                Some(n) => n,
                None => &0,
            };
            let base_n = match cot.get(&'N') {
                Some(n) => n,
                None => &0,
            };

            let total = base_a + base_t + base_g + base_c + base_n;
            let rate_a = *base_a as f64 / total as f64 * 100.0;
            let rate_t = *base_t as f64 / total as f64 * 100.0;
            let rate_g = *base_g as f64 / total as f64 * 100.0;
            let rate_c = *base_c as f64 / total as f64 * 100.0;
            let rate_n = *base_n as f64 / total as f64 * 100.0;

            plot_base_a.insert(pos, rate_a);
            plot_base_t.insert(pos, rate_t);
            plot_base_g.insert(pos, rate_g);
            plot_base_c.insert(pos, rate_c);
            plot_base_n.insert(pos, rate_n);
        }
    }
    let mut dt = vec![];
    dt.push(plot_base_a);
    dt.push(plot_base_t);
    dt.push(plot_base_g);
    dt.push(plot_base_c);
    dt.push(plot_base_n);

    Ok(dt)
}

#[allow(unused_variables)]
#[allow(non_snake_case)]
pub fn open_se(
    r: String,
    out: String,
    figBase: String,
    figQual: String,
    width: u32,
    height: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let fp = File::open(&r)?;
    let reader = BufReader::new(&fp);

    let mut nt_info: hash_base_pos = HashMap::new();
    let mut nt: hash_base = HashMap::new();
    let mut pos_info: hash_qual_pos = HashMap::new();
    let mut pos: hash_qual = HashMap::new();
    let (mut n, mut max_qvalue) = (0u8, 0u8);
    let (mut max_readlen, mut total_seq, total_nt, total_gc) = (0u64, 0u64, 0u64, 0u64);

    for line in reader.lines() {
        n += 1;
        if n == 1 || n == 3 {
            continue;
        }
        if n == 2 {
            if let Ok(seq) = &line {
                total_seq += 1;
                for (order, base) in seq.as_bytes().iter().enumerate() {
                    let real_p = order as u64 + 1;
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
                let length = qual.len() as u64;
                if length > max_readlen {
                    max_readlen = length;
                }

                for (order, qnum) in qual.as_bytes().iter().enumerate() {
                    let q = qnum - 33;
                    if q > max_qvalue {
                        max_qvalue = q;
                    }

                    let real_p = order as u64 + 1;
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
    let df_qual = format_qual(&pos_info, max_readlen, max_qvalue)?;
    let df_base = format_base(&nt_info, max_readlen)?;

    plot_base(df_base.clone(), max_readlen, figBase, width, height)?;
    plot_qual(
        df_qual.clone(),
        max_qvalue,
        max_readlen,
        figQual,
        width,
        height,
    )?;
    show_mat(out, max_qvalue, max_readlen, df_base, df_qual)?;
    Ok(())
}

pub fn show_mat(
    outname: String,
    max_qvalue: u8,
    max_readlen: u64,
    df: Vec<Plotline>,
    qual: Vec<Vec<u64>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut fo = File::create(outname)?;
    let mut head = String::new();
    for m in 0..=max_qvalue {
        if m == 0 {
            head.push_str("Iterm\tA\tT\tG\tC\tN\t");
        }
        if m == max_qvalue {
            head.push_str(format!("{}\n", &m).as_str());
        } else {
            head.push_str(format!("{}\t", &m).as_str());
        }
    }
    fo.write(head.as_bytes())?;

    for i in 0..max_readlen {
        let pos = i + 1;
        let mut nt = format!(
            "pos:{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t",
            pos,
            df[0].get(&pos).unwrap(),
            df[1].get(&pos).unwrap(),
            df[2].get(&pos).unwrap(),
            df[3].get(&pos).unwrap(),
            df[4].get(&pos).unwrap(),
        );
        let mut cot: Vec<String> = vec![];
        for j in qual[i as usize].clone() {
            cot.push(j.to_string());
        }
        nt.push_str(cot.join("\t").as_str());
        nt.push_str("\n");
        fo.write(nt.as_bytes())?;
    }
    Ok(())
}
