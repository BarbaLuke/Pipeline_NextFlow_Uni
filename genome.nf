    canale_genomeDir_FirstPass=Channel.fromPath(params.genomeDir_FirstPass)
    canale_genomeDir_FirstPass2=Channel.fromPath(params.genomeDir_FirstPass)
    canale_genomeDir=Channel.fromPath(params.genomeDir)
    canale_genomeDir2=Channel.fromPath(params.genomeDir)
    canale_runDir_FirstPass=Channel.fromPath(params.runDir_FirstPass)
    canale_runDir=Channel.fromPath(params.runDir)
    canale_genomeFa=Channel.fromPath(params.genomeFa)
    canale_trim=Channel.fromPath(params.fastqDir)


    // --- In questa pipeline andremo a vedere come eseguire 2-pass mapping with re-generated genome con STAR

    // questo è il processo per generare il file di riferimento del genoma
    process generate_genome_indexing_file{

        // in input ho il file in formato .fa di indexing Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
        // e la directory data_dir/reference/genome_index
        input:
        path indexing_file from canale_genomeFa
        path genomeDir_first_pass from canale_genomeDir_FirstPass

        // in output mi serve avere sempre la stessa directory semplicemente per avere una prosecuzione del flow corretto
        output:
        path genomeDir_first_pass_success

        script:
        genomeDir_first_pass_success = genomeDir_first_pass
        
        """
        STAR --runMode genomeGenerate --genomeDir ${genomeDir_first_pass} --genomeFastaFiles ${indexing_file} --runThreadN 50
        """
    }

    // in questo processo eseguiamo il primo passo 
    process first_pass_STAR{

        // per proseguire con il flow corretto metto in input il canale di uscita da generate_genome_indexing_file
        input:
        path genomeDir_first_pass2 from genomeDir_first_pass_success
        val yet_done from done

        // come descritto dai comandi in .sh mandatomi come esempio di pipeline l'unico file in output da questa fase
        // che verrà utilizzato successivamente è SJ.out.tab, lo metto quindi in un canale in output così posso riprenderlo dopo
        output:
        path "SJ.out.tab" into out_firstPass

        // questo piccolo script mi serve per sostituire un comando trovato sempre nel file .sh di esempio di pipeline
        // mi serve avere la lista dei file trimmed (con path relativo) divisi con virgola, tutta su una stringa
        script:
        myDir = file(params.fastqTrimDir)
        file_list = ""
        myDir.eachFile { item ->
        if( item.getExtension() == "fastq.gz" ) {
            if(file_list == ""){
                file_list += params.fastqTrimDir + item.getName()
            }else{
                file_list +="," 
                file_list += params.fastqTrimDir + item.getName()
            }
        }}

        """
        STAR --runMode alignReads --genomeDir ${genomeDir_first_pass2} --readFilesCommand zcat --readFilesIn ${file_list} --runThreadN 50
        """
    }

    // questo è il processo di trimming delle reads di bassa qualità 
    process trimming{
        // questo mi serve per copiare i files di output in una directory a parte nel caso in cui
        // serva averli li dentro invece che dentro la workdir
        publishDir params.fastqTrimDir , mode: 'copy', overwrite: true

        // quando inserisco come canale un filePair devo per forza inserire una tupla in cui il prmo oggetto è un valore
        // il secondo invece sono i files che posso richiamare uno ad uno come in un array (vedere script)
        input:
        path to_trim from canale_trim

        // visto che trimmomatic mi restituisce tanti file, quello che faccio è dichiarare in output  una tupla
        // così da avere a disposizione tutti i file in caso di bisogno
        output:
        path trimmed into canale_file_trimmed
        path trimmed_summary into just_for_saving
        val done
        
        script:
        trimmed = to_trim.getBaseName() + ".trim.fastq"
        trimmed_summary = to_trim.getBaseName() + ".trim.summary"
        // in questo caso richiamo trimmomatic senza bisogno di java perchè ho eseguito un installazione da apt in ubuntu
        """
        TrimmomaticSE -phred33 -threads 3 -quiet ${to_trim} ${trimmed} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -summary ${trimmed_summary}
         """
    }

    // questo è il processo per ri-generare il file di riferimento del genoma
    process REgenerate_genome_indexing_file{

        input:
        path genome_file_fa from canale_genomeFa
        path genomeDir from canale_genomeDir
        path file_out from out_firstPass
        
        output:
        path genomeDir_success

        script:
        genomeDir_success = genomeDir

        """
        STAR --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles ${genome_file_fa} --sjdbFileChrStartEnd ${file_out} --sjdbOverhang 75 --runThreadN 50
        """
    }

    // in questo processo eseguiamo il secondo passo
    process second_pass_STAR{
        publishDir params.runDir , mode: 'copy', overwrite: true

        input:
        path genomeDir2 from genomeDir_success
        path files_trimmed from canale_file_trimmed

        output:
        path file_bam in out_secondPass

        script:
        name = files_trimmed.getSimpleName()
        file_bam = name + "Aligned.sortedByCoord.out.bam"

        """
        STAR --runMode alignReads --genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 --genomeDir ${genomeDir2} --readFilesCommand zcat --readFilesIn ${files_trimmed} --runThreadN 25 --outFileNamePrefix ${name} --outSAMtype BAM SortedByCoordinate
        """
    }
