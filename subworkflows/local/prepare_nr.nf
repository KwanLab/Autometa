
include { DIAMOND_MAKEDB } from './../../modules/local/diamond_makedb.nf'

process DOWNLOAD_NR {
    tag "Downloading nr.gz (>100GB download. May take some time.)"
    label 'process_low'
    label 'process_long'

    println '\033[0;34m Downloading nr.gz from NCBI, this may take a long time. \033[0m'

    conda "conda-forge::rsync=3.2.3"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path("nr.gz")       , emit: singlefile
        path "versions.yml" , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz' 'nr.gz'

        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5' 'nr.gz.md5'

        md5sum -c *.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //' | sed 's/\s.*//')
        END_VERSIONS
        """
}

process TEST_DOWNLOAD {
    // For development work so you don't download the entire nr.gz database
    tag "Downloading small set of FASTA"
    label 'process_low'

    conda "conda-forge::rsync=3.2.3"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path "nr.gz", emit: singlefile

    when:
        task.ext.when == null || task.ext.when

    script:
        """

cat <<-END_VERSIONS > nr
>KJX92028.1 hypothetical protein TI39_contig5958g00003 [Zymoseptoria brevis]
MAWTRQLVPLMLLFCGAHGLQRSSTATDQLSNSALQALGSHADLAAFVNDVEAVPEIANVILAHRGITIMAPVDSAWLRV
DAIKRRNPAFLAWHIMNANVLTSDVPLVQYEQHPGITIPTFLSGSKNWTYSGEPASLISGGQSLTAITLKTEDNVIWVSG
ASNVSYIKQANISYDRGIIHKIDPALQFPTSAYETAFAVGLYSYCWAVFTAGLDQEIRRIPNSTFLLPINEAFHAALPFL
LGASREEFKRIVYRHVIPGRVLWSHEFYNASHETFEGSIVQIRGGNGRRWFVDDAMILDGSDKPLYNGVGHVVNKVLLPT
>EFG1759503.1 decarboxylating NADP(+)-dependent phosphogluconate dehydrogenase [Escherichia coli]EGJ4377881.1 decarboxylating NADP(+)-dependent phosphogluconate dehydrogenase [Escherichia coli]
LKPYLDKGDIIIDGGNTFFQDTIRRNRELSAEGFNFIGTGVSGGEEGALKGPSIMPGGQKEAYELVAPILTKIAAVAEDG
EPCVTYIGADGAGHYVKMVHNGIEYGDMQLIAEAYSLLKGGLNLTNEELAQTFTEWNNGELSSYLIDITKDIFTKKDEDG
NYLVDVILDEAANKGTGKWTSQSALDLGEPLSLITESVFARYISSLKEQRVAASKVLSGPQAQPAGDKGEFIEKVRRALY
LGKIVSYAQGFSQLRAASEEYNWDLNYGEIAKIFRAGCIIRAQFLQKITDAYIENPQIANLLLAPYFKQIADNYQQALRE
VVAYAVQNGIPVPTFAAAVAYYDSYRAAVLPANLIQAQRDYFGAHTYKRIDKEGVFHTEWL
>WP_198835266.1 pilus assembly protein [Paracoccus sp. IB05]MBJ2149627.1 pilus assembly protein [Paracoccus sp. IB05]
MTWRPLQRFLTRSDAAVTAEFVIVFPLVLALIFLIVFISMYISAASDLQQVVHELARYSYRYAGRPEANQLCATLERDAV
PILVNASLLLHPENFTLISCSPPQGPDRIIVITASYDFAGSFVQSVGRTLGLSIGTISRQSLFIP
>MBD3193859.1 hypothetical protein [Candidatus Lokiarchaeota archaeon]MBD3198741.1 hypothetical protein [Candidatus Lokiarchaeota archaeon]
MKKGFIVLILIALVSAGGLILFFYYSNDSGNGNFNTNSEKMIINHNHAHLEDFTSIPSEWIIAAKANLSIVYWHTSHGSQ
ITTGMSLLDAFMGDNDVYEFNNAGTGGALHYHEPSIDYSRRDLTGYTDQFDDETRTFLSSNPEYNVVIWSWCGLDKNNAS
INAYLTNMNQLESEYPNVHFVYMTAHLEGTGEDGDLHIYNQMIRRYCNKNNKTLYDFGDIESYNPENEYFLDRDANDGCY
YDSDGNGSLDANWATEWQSTHDGTHTYPNGGEWYDCSPAHSEAVNGNLKAYAAWYLFARLAGWNGT
>UMM52736.1 protein ORF58 [Lake sturgeon herpesvirus]
MGSMVKKRSRSLIPTSSITRWKTQSLKRPKATCASLRLTPRSTLSPQCHAGYGQSSPGANGLNRPVIDTWTRPSTAFGPS
TSLGWTPQTHIFLNGNFVSHTHGCSPAFFTATQHVNIVYNKKQQTSVFAPHLLPHKQIQSGTVLTDNNKFVTDKKKTFSV
QGVKNTRIEFTSLKNRSSNYTTNCRPLYQPAFQQFFELTGLCHGETSVTMSAMVVNNVNYTTCLYGLTNPFSFNFKICKD
HKKFHNTLFFPSVNLYKQAKGRQHQIFESRYINSQKIYPGDVNQFGFYLQTVVAQTEYDPCLNWYFCRHFEATKSFLNTP
NKTLILWFNERFYLAHPQVDIADPASYWPAYVTFMDLCVTPHLNHFIGFFSSGFGQYHNKNPEFIHLIPFLIFGAARGHN
QGLDLIASYAHRLSRLQRHESLLELRLILQIAVELLKNPQITLCDDPVRGMELSYPQSDDPDNDREKRAKKRRLVVVTKP
LCPPATVVRPLAGHQQSLVKKIQVYCQTCRRG
END_VERSIONS

gzip nr

        """
}

workflow PREPARE_NR_DB {

    main:
        ch_versions = Channel.empty()

        // TODO: this if/else can be simplified
        if (file("${params.nr_dmnd_dir}/nr.dmnd").exists()){
            // skip huge download and db creation if nr.dmnd already exists
            out_ch = file("${params.nr_dmnd_dir}/nr.dmnd")
        } else if (file("${params.nr_dmnd_dir}/nr.gz").exists()){
            // skip huge download if nr.gz already exists
            DIAMOND_MAKEDB(file("${params.nr_dmnd_dir}/nr.gz"), "nr")
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
            out_ch = DIAMOND_MAKEDB.out.diamond_db

        } else if (params.debug){
            TEST_DOWNLOAD().singlefile
                .set{nr_db_ch}

            DIAMOND_MAKEDB(nr_db_ch, "nr")
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
            out_ch = DIAMOND_MAKEDB.out.diamond_db

        } else if (params.large_downloads_permission) {
            DOWNLOAD_NR().singlefile
                .set{nr_db_ch}
            ch_versions = ch_versions.mix(DOWNLOAD_NR.out.versions)
            DIAMOND_MAKEDB(nr_db_ch, "nr")
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
            out_ch = DIAMOND_MAKEDB.out.diamond_db

        } else {
            println '\033[0;34m Neither nr.dmnd or nr.gz were found and `--large_downloads_permission` is set to false. \033[0m'
        }

    emit:
        diamond_db = out_ch
        versions = ch_versions
}
