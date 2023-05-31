use std::fs::File;
use std::io::{ BufReader, BufRead};


pub struct GenomicInterval {
    chromosome: String,
    start: u32,
    end: u32,
}

pub struct ChromSizes {
    chr1: u32,
    chr2: u32,
    chr3: u32,
    chr4: u32,
    chr5: u32,
    chr6: u32,
    chr7: u32,
    chr8: u32,
    chr9: u32,
    chr10: u32,
    chr11: u32,
    chr12: u32,
    chr13: u32,
    chr14: u32,
    chr15: u32,
    chr16: u32,
    chr17: u32,
    chr18: u32,
    chr19: u32,
    chr20: u32,
    chr21: u32,
    chr22: u32,
    chr_x: u32,
    chr_y: u32,
    chr_m: u32,
}

pub fn parse_bed(file_path: &str) -> Vec<GenomicInterval> {

    let file = File::open(file_path).expect("Failed to open BED file");
    let reader = BufReader::new(file);

    let mut intervals = Vec::new();

    for line in reader.lines() {
        if let Ok(line) = line {
            let fields: Vec<&str> = line.split('\t').collect();

            let chromosome = fields[0].to_string();
            let start = fields[1].parse::<u32>().expect("Failed to parse start position");
            let end = fields[2].parse::<u32>().expect("Failed to parse end position");

            let interval = GenomicInterval {
                chromosome,
                start,
                end,
            };

            intervals.push(interval);
        }
    }
    intervals
}

pub fn read_chrom_sizes(file_path: &str) -> ChromSizes {
    let file = File::open(file_path).expect("Failed to open chrom.sizes file");
    let reader = BufReader::new(file);

    let mut chrom_sizes = ChromSizes {
        chr1: 0,
        chr2: 0,
        chr3: 0,
        chr4: 0,
        chr5: 0,
        chr6: 0,
        chr7: 0,
        chr8: 0,
        chr9: 0,
        chr10: 0,
        chr11: 0,
        chr12: 0,
        chr13: 0,
        chr14: 0,
        chr15: 0,
        chr16: 0,
        chr17: 0,
        chr18: 0,
        chr19: 0,
        chr20: 0,
        chr21: 0,
        chr22: 0,
        chr_x: 0,
        chr_y: 0,
        chr_m: 0,
    };

    for line in reader.lines() {
        if let Ok(line) = line {
            let fields: Vec<&str> = line.split('\t').collect();

            let chromosome = fields[0].to_string();
            let size = fields[1].parse::<u32>().expect("Failed to parse chromosome size");

            match chromosome.as_str() {
                "chr1" => chrom_sizes.chr1 = size,
                "chr2" => chrom_sizes.chr2 = size,
                "chr3" => chrom_sizes.chr3 = size,
                "chr4" => chrom_sizes.chr4 = size,
                "chr5" => chrom_sizes.chr5 = size,
                "chr6" => chrom_sizes.chr6 = size,
                "chr7" => chrom_sizes.chr7 = size,
                "chr8" => chrom_sizes.chr8 = size,
                "chr9" => chrom_sizes.chr9 = size,
                "chr10" => chrom_sizes.chr10 = size,
                "chr11" => chrom_sizes.chr11 = size,
                "chr12" => chrom_sizes.chr12 = size,
                "chr13" => chrom_sizes.chr13 = size,
                "chr14" => chrom_sizes.chr14 = size,
                "chr15" => chrom_sizes.chr15 = size,
                "chr16" => chrom_sizes.chr16 = size,
                "chr17" => chrom_sizes.chr17 = size,
                "chr18" => chrom_sizes.chr18 = size,
                "chr19" => chrom_sizes.chr19 = size,
                "chr20" => chrom_sizes.chr20 = size,
                "chr21" => chrom_sizes.chr21 = size,
                "chr22" => chrom_sizes.chr22 = size,
                "chrX" => chrom_sizes.chr_x = size,
                "chrY" => chrom_sizes.chr_y = size,
                "chrM" => chrom_sizes.chr_m = size,
                _ => println!("Ignoring unknown chromosome: {}", chromosome),
            }
        }   
    }

    // check all values greater than 0
    assert!(chrom_sizes.chr1 > 0);
    assert!(chrom_sizes.chr2 > 0);
    assert!(chrom_sizes.chr3 > 0);
    assert!(chrom_sizes.chr4 > 0);
    assert!(chrom_sizes.chr5 > 0);
    assert!(chrom_sizes.chr6 > 0);
    assert!(chrom_sizes.chr7 > 0);
    assert!(chrom_sizes.chr8 > 0);
    assert!(chrom_sizes.chr9 > 0);
    assert!(chrom_sizes.chr10 > 0);
    assert!(chrom_sizes.chr11 > 0);
    assert!(chrom_sizes.chr12 > 0);
    assert!(chrom_sizes.chr13 > 0);
    assert!(chrom_sizes.chr14 > 0);
    assert!(chrom_sizes.chr15 > 0);
    assert!(chrom_sizes.chr16 > 0);
    assert!(chrom_sizes.chr17 > 0);
    assert!(chrom_sizes.chr18 > 0);
    assert!(chrom_sizes.chr19 > 0);
    assert!(chrom_sizes.chr20 > 0);
    assert!(chrom_sizes.chr21 > 0);
    assert!(chrom_sizes.chr22 > 0);
    assert!(chrom_sizes.chr_x > 0);
    assert!(chrom_sizes.chr_y > 0);
    assert!(chrom_sizes.chr_m > 0);

    chrom_sizes

}

/// line up all chromosomes in a single vector, then convert to binary using intervals
pub fn make_binary_vector(intervals: &[GenomicInterval], chrom_sizes: &ChromSizes) -> Vec<u8> {
    let mut binary_vector = Vec::new();

    // for chrom_size in chrom_sizes, add 0s to binary_vector
    let mut chr1_binary_vector = vec![0; chrom_sizes.chr1 as usize];
    let mut chr2_binary_vector = vec![0; chrom_sizes.chr2 as usize];
    let mut chr3_binary_vector = vec![0; chrom_sizes.chr3 as usize];
    let mut chr4_binary_vector = vec![0; chrom_sizes.chr4 as usize];
    let mut chr5_binary_vector = vec![0; chrom_sizes.chr5 as usize];
    let mut chr6_binary_vector = vec![0; chrom_sizes.chr6 as usize];
    let mut chr7_binary_vector = vec![0; chrom_sizes.chr7 as usize];
    let mut chr8_binary_vector = vec![0; chrom_sizes.chr8 as usize];
    let mut chr9_binary_vector = vec![0; chrom_sizes.chr9 as usize];
    let mut chr10_binary_vector = vec![0; chrom_sizes.chr10 as usize];
    let mut chr11_binary_vector = vec![0; chrom_sizes.chr11 as usize];
    let mut chr12_binary_vector = vec![0; chrom_sizes.chr12 as usize];
    let mut chr13_binary_vector = vec![0; chrom_sizes.chr13 as usize];
    let mut chr14_binary_vector = vec![0; chrom_sizes.chr14 as usize];
    let mut chr15_binary_vector = vec![0; chrom_sizes.chr15 as usize];
    let mut chr16_binary_vector = vec![0; chrom_sizes.chr16 as usize];
    let mut chr17_binary_vector = vec![0; chrom_sizes.chr17 as usize];
    let mut chr18_binary_vector = vec![0; chrom_sizes.chr18 as usize];
    let mut chr19_binary_vector = vec![0; chrom_sizes.chr19 as usize];
    let mut chr20_binary_vector = vec![0; chrom_sizes.chr20 as usize];
    let mut chr21_binary_vector = vec![0; chrom_sizes.chr21 as usize];
    let mut chr22_binary_vector = vec![0; chrom_sizes.chr22 as usize];
    let mut chr_x_binary_vector = vec![0; chrom_sizes.chr_x as usize];
    let mut chr_y_binary_vector = vec![0; chrom_sizes.chr_y as usize];
    let mut chr_m_binary_vector = vec![0; chrom_sizes.chr_m as usize];

    // for each interval, add 1s to binary_vector
    for interval in intervals {
        match interval.chromosome.as_str() {
            "chr1" => {
                for i in interval.start..interval.end {
                    chr1_binary_vector[i as usize] = 1;
                }
            },
            "chr2" => {
                for i in interval.start..interval.end {
                    chr2_binary_vector[i as usize] = 1;
                }
            },
            "chr3" => {
                for i in interval.start..interval.end {
                    chr3_binary_vector[i as usize] = 1;
                }
            },
            "chr4" => {
                for i in interval.start..interval.end {
                    chr4_binary_vector[i as usize] = 1;
                }
            },
            "chr5" => {
                for i in interval.start..interval.end {
                    chr5_binary_vector[i as usize] = 1;
                }
            },
            "chr6" => {
                for i in interval.start..interval.end {
                    chr6_binary_vector[i as usize] = 1;
                }
            },
            "chr7" => {
                for i in interval.start..interval.end {
                    chr7_binary_vector[i as usize] = 1;
                }
            },
            "chr8" => {
                for i in interval.start..interval.end {
                    chr8_binary_vector[i as usize] = 1;
                }
            },
            "chr9" => {
                for i in interval.start..interval.end {
                    chr9_binary_vector[i as usize] = 1;
                }
            },
            "chr10" => {
                for i in interval.start..interval.end {
                    chr10_binary_vector[i as usize] = 1;
                }
            },
            "chr11" => {
                for i in interval.start..interval.end {
                    chr11_binary_vector[i as usize] = 1;
                }
            },
            "chr12" => {
                for i in interval.start..interval.end {
                    chr12_binary_vector[i as usize] = 1;
                }
            },
            "chr13" => {
                for i in interval.start..interval.end {
                    chr13_binary_vector[i as usize] = 1;
                }
            },
            "chr14" => {
                for i in interval.start..interval.end {
                    chr14_binary_vector[i as usize] = 1;
                }
            },
            "chr15" => {
                for i in interval.start..interval.end {
                    chr15_binary_vector[i as usize] = 1;
                }
            },
            "chr16" => {
                for i in interval.start..interval.end {
                    chr16_binary_vector[i as usize] = 1;
                }
            },
            "chr17" => {
                for i in interval.start..interval.end {
                    chr17_binary_vector[i as usize] = 1;
                }
            },
            "chr18" => {
                for i in interval.start..interval.end {
                    chr18_binary_vector[i as usize] = 1;
                }
            },
            "chr19" => {
                for i in interval.start..interval.end {
                    chr19_binary_vector[i as usize] = 1;
                }
            },
            "chr20" => {
                for i in interval.start..interval.end {
                    chr20_binary_vector[i as usize] = 1;
                }
            },
            "chr21" => {
                for i in interval.start..interval.end {
                    chr21_binary_vector[i as usize] = 1;
                }
            },
            "chr22" => {
                for i in interval.start..interval.end {
                    chr22_binary_vector[i as usize] = 1;
                }
            },
            "chrX" => {
                for i in interval.start..interval.end {
                    chr_x_binary_vector[i as usize] = 1;
                }
            },
            "chrY" => {
                for i in interval.start..interval.end {
                    chr_y_binary_vector[i as usize] = 1;
                }
            },
            "chrM" => {
                for i in interval.start..interval.end {
                    chr_m_binary_vector[i as usize] = 1;
                }
            },
            _ => println!("Ignoring unknown chromosome: {}", interval.chromosome),
        }
    }

    // add all chromosomes to binary_vector
    binary_vector.append(&mut chr1_binary_vector);
    binary_vector.append(&mut chr2_binary_vector);
    binary_vector.append(&mut chr3_binary_vector);
    binary_vector.append(&mut chr4_binary_vector);
    binary_vector.append(&mut chr5_binary_vector);
    binary_vector.append(&mut chr6_binary_vector);
    binary_vector.append(&mut chr7_binary_vector);
    binary_vector.append(&mut chr8_binary_vector);
    binary_vector.append(&mut chr9_binary_vector);
    binary_vector.append(&mut chr10_binary_vector);
    binary_vector.append(&mut chr11_binary_vector);
    binary_vector.append(&mut chr12_binary_vector);
    binary_vector.append(&mut chr13_binary_vector);
    binary_vector.append(&mut chr14_binary_vector);
    binary_vector.append(&mut chr15_binary_vector);
    binary_vector.append(&mut chr16_binary_vector);
    binary_vector.append(&mut chr17_binary_vector);
    binary_vector.append(&mut chr18_binary_vector);
    binary_vector.append(&mut chr19_binary_vector);
    binary_vector.append(&mut chr20_binary_vector);
    binary_vector.append(&mut chr21_binary_vector);
    binary_vector.append(&mut chr22_binary_vector);
    binary_vector.append(&mut chr_x_binary_vector);
    binary_vector.append(&mut chr_y_binary_vector);
    binary_vector.append(&mut chr_m_binary_vector);

    binary_vector
    
}

fn compute_dice_coeff(binary_vector1: &[u8], binary_vector2: &[u8]) -> f64 {

    // first ensure nothing weird happend
    assert_eq!(binary_vector1.len(), binary_vector2.len());

    let mut intersection = 0;
    let mut union = 0;

    for i in 0..binary_vector1.len() {
        if binary_vector1[i] == 1 && binary_vector2[i] == 1 {
            intersection += 1;
        }
        if binary_vector1[i] == 1 || binary_vector2[i] == 1 {
            union += 1;
        }
    }

    let dice_coeff = (2.0 * intersection as f64) / (union as f64);

    dice_coeff
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed() {
        let bed = "test_data/test.bed";
        let intervals = parse_bed(bed);
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].chromosome, "chr1");
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
        assert_eq!(intervals[1].chromosome, "chr2");
        assert_eq!(intervals[1].start, 300);
        assert_eq!(intervals[1].end, 400);
    }

    #[test]
    fn test_read_chrom_sizes() {
        let sizes_file = "test_data/hg38.chrom.sizes";
        let chrom_sizes = read_chrom_sizes(sizes_file);

        assert_eq!(chrom_sizes.chr1, 248956422);
    }
    #[test]
    fn test_create_binary_vector() {
        let sizes_file = "test_data/hg38.chrom.sizes";
        let chrom_sizes = read_chrom_sizes(sizes_file);

        let bed = "test_data/test.bed";
        let intervals = parse_bed(bed);

        let binary_vector = make_binary_vector(&intervals, &chrom_sizes);
        assert!(binary_vector.len() > 0)
    }
    #[test]
    fn test_coeff_full() {
        let sizes_file = "test_data/hg38.chrom.sizes";
        let chrom_sizes = read_chrom_sizes(sizes_file);


        let bed = "test_data/test.bed";
        let intervals = parse_bed(bed);

        let binary_vector1 = make_binary_vector(&intervals, &chrom_sizes);
        let binary_vector2 = make_binary_vector(&intervals, &chrom_sizes);

        let dice_coeff = compute_dice_coeff(&binary_vector1, &binary_vector2);
        assert_eq!(dice_coeff, 2.0);
    }
}
