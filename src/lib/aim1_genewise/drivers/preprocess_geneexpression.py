def preprocess():
    pass
    # Decompress gtex data

    # Load as SummarizedExperiment
    # from https://rdrr.io/bioc/slinky/man/coerce.html
    # for build/demo only.  You MUST use your own key when using the slinky
    # package.
    user_key < - httr::content(httr::GET('https://api.clue.io/temp_api_key'),
    as='parsed')$user_key
    sl < - Slinky(user_key,
                  system.file('extdata', 'demo.gctx',
                              package='slinky'),
                  system.file('extdata', 'demo_inst_info.txt',
                              package='slinky'))
    sumex < - as(sl[, 1:20], "SummarizedExperiment")

    # Filter
