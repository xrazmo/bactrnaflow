nextflow.enable.dsl=2

params.meta_csv=""
params.ref_dir=""
params.fq_dir=""
params.output=""
params.contig_dir=""
params.prodigal_tf="$baseDir/assets/Escherichia_coli.trn"
params.dmnd_dir = "$baseDir/assets/dmnd"

include {FASTQC} from "$baseDir/modules/nf-core/fastqc/main"
include {TRIMMOMATIC} from "$baseDir/modules/nf-core/trimmomatic/main"
include {STAR_GENOMEGENERATE} from "$baseDir/modules/nf-core/star/genomegenerate/main"
include {STAR_ALIGN} from "$baseDir/modules/nf-core/star/align/main"
include {SUBREAD_FEATURECOUNTS} from "$baseDir/modules/nf-core/subread/featurecounts/main"

include {PROKKA} from "$baseDir/modules/local/prokka/main"
// include {ROARY} from "$baseDir/modules/local/roary/main"
include {DIAMOND_BLASTX} from "$baseDir/modules/local/diamond/blastx/main"

workflow PrepRef{

    ref_contigs_ch = Channel.fromPath("${params.contig_dir}/*.gz").map{it->[[id:it.simpleName],it]}
   
    PROKKA(ref_contigs_ch,[],params.prodigal_tf,"${params.output}/prokka")
   
    //ROARY(PROKKA.out.gff.map{it->it[1]}.flatten().collect(),"${params.output}/roary")
    // orf_ch = Channel.fromPath("${params.output}/prokka/*/*.ffn").map{it->[[id:it.simpleName],it]}
   

    dmnd_ch = Channel.fromPath("${params.dmnd_dir}/*.dmnd")
    def diamond_cols = "qseqid sseqid pident qcovhsp scovhsp mismatch gaps evalue bitscore length qlen slen qstart qend sstart send stitle"
    
    PROKKA.out.ffn.combine(dmnd_ch).combine(Channel.from('txt')).combine(Channel.from(tuple(diamond_cols)))
                    .map{it -> [[id:it[0].id+"__"+it[2].simpleName],it[1],it[2],it[3],it[4]]}
                    .multiMap {it -> fasta: tuple(it[0],it[1])
                                    db: it[2]
                                    ext: it[3]
                                    col: it[4]}
                    .set{diaParam}
    
    
    DIAMOND_BLASTX(diaParam.fasta,diaParam.db,diaParam.ext,diaParam.col)
    DIAMOND_BLASTX.out.txt.map{it-> it[1]}.flatten().collectFile(storeDir:"${params.output}/prokka/extra_annotation")                
}


workflow {

    fq_ch = Channel.fromFilePairs("${params.fq_dir}/*_{R1,R2}.fastq.gz")
                    .map{it->[[id:it[0],single_end:false],it[1]]}

    ref_ch = Channel.fromFilePairs("${params.ref_dir}/*/*.{fna,gtf}")
                    .multiMap{it-> 
                                genomes:[[id:it[0]],it[1]]
                                gtfs: [[id:it[0]],it[1][1]]
                                fnas:[[id:it[0]],it[1][0]]
                            }

    // CSV file linking Fastq files with their corresponding reference genomes    
    meta=Channel.fromPath(params.meta_csv)
                .splitCsv(header:true)
                .map{it->[[id:it.fq_id,single_end:false],[id:it.ref_id]]}                

    FASTQC(fq_ch)

    TRIMMOMATIC(fq_ch)

    STAR_GENOMEGENERATE(ref_ch.genomes)
    
    // creating a tuple for STAR_ALIGN [ids,PE-reads,genome index, gtf file]
     
     ch1 = TRIMMOMATIC.out.trimmed_reads.join(meta).map{it->[it[2],it[0],it[1]]}
     ch2 = ref_ch.gtfs.join(STAR_GENOMEGENERATE.out.index)

     star_input_ch = ch2.cross(ch1).map{it->[it[1][1]+['ref_id':it[0][0].id],it[1][2],it[0][2],it[0][1]]}
        
   
    STAR_ALIGN(star_input_ch,false,'','')


    fcount_input_ch = ref_ch.gtfs.cross(STAR_ALIGN.out.bam
                 .map{it->[[id:it[0].ref_id],it[0],it[1]]}).map{it->[it[1][1],it[1][2],it[0][1]]}

     
    
    SUBREAD_FEATURECOUNTS(fcount_input_ch)

    FASTQC.out.html.map{it->it[1]}.flatten().collectFile(storeDir:"${params.output}/fastqc")
    SUBREAD_FEATURECOUNTS.out.counts.map{it->it[1]}.flatten().collectFile(storeDir:"${params.output}/feature_count")
    SUBREAD_FEATURECOUNTS.out.summary.map{it->it[1]}.flatten().collectFile(storeDir:"${params.output}/feature_count")
}