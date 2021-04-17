canale_trim=Channel.fromPath(params.fastqDir)
canale_fastQC_before=Channel.fromPath(params.fastqDir)
canale_whereQC_before=Channel.fromPath(params.fastQC_before)
canale_whereQC_after=Channel.fromPath(params.fastQC_after)

// questo è il processo del quality-control prima del trimming
process fastqc_before{
    // per questo processo l'output viene direttamente inserito in una cartella
    // questo poichè viene generata una pagina HTML da visionare 
    // tutti il necessario è nella cartella QC/before dentro data_dir

    input:
    path to_control_before from canale_fastQC_before
    path where_is_before from canale_whereQC_before

    """
    fastqc -o ${where_is_before} ${to_control_before}
    """
}

// questo è il processo di trimming delle reads di bassa qualità 
process trimming{
    // questo mi serve per copiare i files di output in una directory per a parte nel caso in cui
    // serva averli li dentro invece che dentro la workdir
    publishDir params.fastqTrimDir , mode: 'copy', overwrite: true

    // quando inserisco come canale un filePair devo per forza inserire una tupla in cui il prmo oggetto è un valore
    // il secondo invece sono i files che posso richiamare uno ad uno come in un array (vedere script)
    input:
    path to_trim from canale_trim

    // visto che trimmomatic mi restituisce tanti file, quello che faccio è dichiarare in output  una tupla
    // così da avere a disposizione tutti i file in caso di bisogno
    output:
    path trimmed into canale_fastQC_after
    path trimmed_summary into solo_per_salvarlo

    script:
    trimmed = to_trim.getBaseName() + "_trimmed.fastq"
    trimmed_summary = to_trim.getBaseName() + ".trim.summary"

    // in questo caso richiamo trimmomatic senza bisogno di java perchè ho eseguito un installazione da apt in ubuntu
    """
    TrimmomaticSE -phred33 -threads 3 -quiet ${to_trim} ${trimmed} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -summary ${trimmed_summary}
    """
}

// questo è nuovamente il processo di quality-control però sul file fastaq filtrato delle basi di bassa qualita'
process fastqc_after{
    // ripeto il quality control per vedere cosa è cambiato con il trimming
    // il processo è del tutto uguale a prima se non fosse che prendi i file dal canale creato dall'output del processo di trimming
    // tutti il necessario è nella cartella QC/after dentro data_dir

    input:
    path to_control_after from canale_fastQC_after
    path where_is_after from canale_whereQC_after

    """
    fastqc -o ${where_is_after} ${to_control_after}
    """
}
