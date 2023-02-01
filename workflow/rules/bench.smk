from pathlib import Path


out = "results/bench/"

rule pgen_write:
    params:
        variants = config["num_variants"],
        samples = config["num_samples"],
    output:
        pgen = temp(out+"write/pgen/{rep}.pgen"),
        pvar = temp(out+"write/pgen/{rep}.pvar"),
        psam = temp(out+"write/pgen/{rep}.psam"),
        time = out+"write/pgen/{rep}.txt"
    resources:
        runtime="0:20:00",
        queue="hotel",
    log:
        out+"log/write/pgen/{rep}.log"
    benchmark:
        out+"bench/write/pgen/{rep}.txt"
    threads: 5
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/write_fake_gts.py --plink -v {params.variants} "
        "-s {params.samples} {output.pgen} > {output.time} 2> {log}"

rule vcf_write:
    params:
        variants = config["num_variants"],
        samples = config["num_samples"],
    output:
        vcf = temp(out+"write/vcf/{rep}.vcf"),
        time = out+"write/vcf/{rep}.txt"
    resources:
        runtime="10:00:00",
        queue="hotel",
    log:
        out+"log/write/vcf/{rep}.log"
    benchmark:
        out+"bench/write/vcf/{rep}.txt"
    threads: 4
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/write_fake_gts.py -v {params.variants} "
        "-s {params.samples} {output.vcf} > {output.time} 2> {log}"

rule bcf_write:
    params:
        variants = config["num_variants"],
        samples = config["num_samples"],
    output:
        vcf = temp(out+"write/bcf/{rep}.bcf"),
        time = out+"write/bcf/{rep}.txt"
    resources:
        runtime="10:00:00",
        queue="hotel",
    log:
        out+"log/write/bcf/{rep}.log"
    benchmark:
        out+"bench/write/bcf/{rep}.txt"
    threads: 4
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/write_fake_gts.py -v {params.variants} "
        "-s {params.samples} {output.vcf} > {output.time} 2> {log}"

rule pgen_read:
    input:
        pgen = rules.pgen_write.output.pgen,
        pvar = rules.pgen_write.output.pvar,
        psam = rules.pgen_write.output.psam,
    output:
        time = out+"read/pgen/{rep}.txt"
    resources:
        runtime="0:20:00",
        queue="hotel",
    log:
        out+"log/read/pgen/{rep}.log"
    benchmark:
        out+"bench/read/pgen/{rep}.txt"
    threads: 5
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/read_gts.py {input.pgen} > {output.time} 2> {log}"

rule vcf_read:
    input:
        vcf = rules.vcf_write.output.vcf,
    output:
        time = out+"read/vcf/{rep}.txt"
    resources:
        runtime="10:00:00",
        queue="hotel",
    log:
        out+"log/read/vcf/{rep}.log"
    benchmark:
        out+"bench/read/vcf/{rep}.txt"
    threads: 4
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/read_gts.py {input.vcf} > {output.time} 2> {log}"

rule bcf_read:
    input:
        vcf = rules.bcf_write.output.vcf,
    output:
        time = out+"read/bcf/{rep}.txt"
    resources:
        runtime="10:00:00",
        queue="hotel",
    log:
        out+"log/read/bcf/{rep}.log"
    benchmark:
        out+"bench/read/bcf/{rep}.txt"
    threads: 4
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/read_gts.py {input.vcf} > {output.time} 2> {log}"

rule report:
    input:
        write_pgen = expand(rules.pgen_write.output.time, rep=range(config["reruns"])),
        write_vcf = expand(rules.vcf_write.output.time, rep=range(config["reruns"])),
        write_bcf = expand(rules.bcf_write.output.time, rep=range(config["reruns"])),
        read_pgen = expand(rules.pgen_read.output.time, rep=range(config["reruns"])),
        read_vcf = expand(rules.vcf_read.output.time, rep=range(config["reruns"])),
        read_bcf = expand(rules.bcf_read.output.time, rep=range(config["reruns"])),
    output:
        report = out+"report.txt"
    resources:
        runtime="0:01:00"
    log:
        out+"log/report.log"
    benchmark:
        out+"bench/report.txt"
    threads: 4
    conda:
        "../envs/haptools.yml"
    shell:
        "{{ echo -e 'write_pgen\\tread_pgen\\twrite_vcf\\tread_vcf\\twrite_bcf\\tread_bcf'; paste "
        "<(cat {input.write_pgen}) <(cat {input.read_pgen})"
        " <(cat {input.write_vcf}) <(cat {input.read_vcf}) "
        " <(cat {input.write_bcf}) <(cat {input.read_bcf})"
        "; }} > {output.report} 2> {log}"
