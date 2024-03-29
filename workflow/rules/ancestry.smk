from pathlib import Path


out = "results/ancestry/"
# get the ID of the causal SNP from the hap filename
config["hap"] = Path(config["hap"])
snp_id = config["hap"].with_suffix("")
if snp_id.suffix == ".hap":
    snp_id = snp_id.with_suffix("")
snp_id = str(snp_id.name)
# create a region param that encodes a 1 Mbp region around the given site
config["region"] = tuple(map(int, config["snp"].split(":")))
config["region"] = tuple(map(str, (
    config["region"][0],
    int(config["region"][1] - 1e6),
    int(config["region"][1] + 1e6),
)))


rule sim_gts:
    input:
        ref = config["reference"],
        samps = config["sample_info"],
        model = lambda wildcards: config["models"][wildcards.samp],
        mapdir = config["mapdir"],
    params:
        chroms = config["region"][0],
        region = config["region"][0]+":"+"-".join(config["region"][1:]),
        out_prefix = lambda w, output: Path(output.gts).with_suffix(""),
        nsamps = 10*10000,
    output:
        gts = temp(out+"sim_gts/{samp}.vcf"),
        bkpt = out+"sim_gts/{samp}.bp",
    resources:
        runtime="3:00:00"
    log:
        out+"logs/sim_gts/{samp}.log"
    benchmark:
        out+"bench/sim_gts/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools simgenotype --invcf {input.ref} --sample_info {input.samps} "
        "--model {input.model} --mapdir {input.mapdir} --chroms {params.chroms} "
        "--out {params.out_prefix} --region {params.region} --popsize {params.nsamps}"
        "&> {log}"

rule vcf2pgen:
    input:
        gts = rules.sim_gts.output.gts,
    params:
        prefix = lambda w, output: Path(output.pgen).with_suffix(""),
    output:
        log = temp(out+"sim_gts/{samp}.log"),
        pgen = out+"sim_gts/{samp}.pgen",
        pvar = out+"sim_gts/{samp}.pvar",
        psam = out+"sim_gts/{samp}.psam",
    resources:
        runtime="0:05:00"
    log:
        out+"logs/vcf2pgen/{samp}.log"
    benchmark:
        out+"bench/vcf2pgen/{samp}.txt"
    threads: 12
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --vcf {input.gts} --make-pgen --snps-only --max-alleles 2 "
        "--threads {threads} --out {params.prefix} &> {log}"

rule transform:
    input:
        pgen = rules.vcf2pgen.output.pgen,
        pvar = rules.vcf2pgen.output.pvar,
        psam = rules.vcf2pgen.output.psam,
        bkpt = rules.sim_gts.output.bkpt,
        hap = config["hap"],
    params:
        region = config["region"][0]+":"+"-".join(config["region"][1:]),
    output:
        pgen = temp(out+"transform/{samp}.pgen"),
        pvar = temp(out+"transform/{samp}.pvar"),
        psam = temp(out+"transform/{samp}.psam"),
    resources:
        runtime="0:04:00"
    log:
        out+"logs/transform/{samp}.log"
    benchmark:
        out+"bench/transform/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools transform -v INFO --ancestry -o {output.pgen} "
        "--region {params.region} {input.pgen} {input.hap} &> {log}"

def sim_pts_input(wildcards):
    source = rules.transform.input
    if wildcards.type == "ancestry":
        source = rules.transform.output
    return {
        "pgen": source.pgen,
        "pvar": source.pvar,
        "psam": source.psam,
    }

rule sim_pts:
    input:
        unpack(sim_pts_input),
        hap = config["hap"],
    params:
        beta = lambda wildcards: wildcards.beta,
        p = lambda w: ["", f"-p {config['prevalence']}"][w.cc == "cc"],
        reps = config["replicates"],
    output:
        pts = out+"sim_pts/{cc}/b{beta}/{type}/{samp}.pheno",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/sim_pts/{cc}/b{beta}/{type}/{samp}.log"
    benchmark:
        out+"bench/sim_pts/{cc}/b{beta}/{type}/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools simphenotype -r {params.reps} -o {output.pts} -v DEBUG {input.pgen} "
        "{params.p} <( zcat {input.hap} | sed 's/\\t0.99$/\\t{params.beta}/' ) &>{log}"

rule gwas:
    input:
        pgen = rules.transform.input.pgen,
        pvar = rules.transform.input.pvar,
        psam = rules.transform.input.psam,
        pts = rules.sim_pts.output.pts,
    params:
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
        cc = lambda wildcards: "--1 " if wildcards.cc == "cc" else ""
    output:
        log = temp(out+"sim_pts/{cc}/b{beta}/{type}/{samp}/{samp}.log"),
        linear = directory(out+"sim_pts/{cc}/b{beta}/{type}/{samp}"),
    resources:
        runtime="0:04:00"
    log:
        out+"logs/gwas/{cc}/b{beta}/{type}/{samp}.log"
    benchmark:
        out+"bench/gwas/{cc}/b{beta}/{type}/{samp}.txt"
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --linear allow-no-covars --variance-standardize {params.cc}"
        "--pheno {input.pts} --pfile {params.in_prefix} --out {params.out_prefix} "
        "--threads {threads} &>{log}"

rule manhattan:
    input:
        linear = expand(
            rules.gwas.output.linear,
            type="ancestry",
            samp=config["models"].keys(),
            allow_missing=True,
        ),
    params:
        linear = lambda wildcards, input: [
            f"-l {i}/{samp}."+snp_id+".glm.linear"
            for i, samp in zip(input.linear, config["models"].keys())
        ],
        snp = snp_id,
        size = "small",
    output:
        png = out+"sim_pts/{cc}/b{beta}/manhattan.pdf",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/manhattan/{cc}/b{beta}/manhattan.log"
    benchmark:
        out+"bench/manhattan/{cc}/b{beta}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py --{params.size} -o {output.png} {params.linear}"
        " -i {params.snp} --no-label &>{log}"

ancestry_pat = expand(
    rules.gwas.output.linear,
    type = "ancestry",
    samp = "Admixed",
    allow_missing=True,
)[0]
normal_pat = expand(
    rules.gwas.output.linear,
    type = "normal",
    samp = "Admixed",
    allow_missing=True,
)[0]

rule power:
    input:
        ancestry = expand(ancestry_pat, beta = config["betas"], allow_missing=True),
        normal = expand(normal_pat, beta = config["betas"], allow_missing=True),
    params:
        binary = lambda wildcards: "--binary " if wildcards.cc == "cc" else "",
        ancestry = lambda w: expand(ancestry_pat, cc=w.cc, allow_missing=True),
        normal = lambda w: expand(normal_pat, cc=w.cc, allow_missing=True),
        snp = snp_id,
    output:
        png = out+"sim_pts/{cc}/power.pdf"
    resources:
        runtime="0:05:00"
    log:
        out+"logs/power/{cc}/log.log"
    benchmark:
        out+"bench/power/{cc}/bench.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/power.py -o {output.png} {params.binary} {params.snp} "
        "'{params.ancestry}' '{params.normal}' &>{log}"
