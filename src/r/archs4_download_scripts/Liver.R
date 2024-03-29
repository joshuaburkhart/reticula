# R script to download selected samples
# Copy code and run on a local machine to initiate download

# Check for dependencies and install if missing
packages <- c("rhdf5")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    print("Install required packages")
    source("https://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
}
library("rhdf5")
library("tools")

destination_file = "human_matrix_download.h5"
extracted_expression_file = "Liver_expression_matrix.tsv"
url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"

# Check if gene expression file was already downloaded and check integrity, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    download.file(url, destination_file, quiet = FALSE)
} else{
    print("Verifying file integrity...")
    checksum = md5sum(destination_file)
    
    if(destination_file == "human_matrix_download.h5"){
        # human checksum (checksum is for latest version of ARCHS4 data)
        correct_checksum = "34197866d7841cc4fb31e09195faa150"
    } else{
        # mouse checksum (checksum is for latest version of ARCHS4 data)
        correct_checksum = "55441d1af9da82c6f3d368c8fa554d42"
    }
    
    if(checksum != correct_checksum){
        print("Existing file looks corrupted or is out of date. Downloading compressed gene expression matrix again.")
        download.file(url, destination_file, quiet = FALSE)
    } else{
        print("Latest ARCHS4 file already exists.")
    }
}

checksum = md5sum(destination_file)
if(destination_file == "human_matrix_download.h5"){
    # human checksum (checksum is for latest version of ARCHS4 data)
    correct_checksum = "34197866d7841cc4fb31e09195faa150"
} else{
    # mouse checksum (checksum is for latest version of ARCHS4 data)
    correct_checksum = "55441d1af9da82c6f3d368c8fa554d42"
}

if(checksum != correct_checksum){
    print("File download ran into problems. Please try to download again. The files are also available for manual download at http://amp.pharm.mssm.edu/archs4/download.html.")
} else{
    # Selected samples to be extracted
    samp = c("GSM742944","GSM2326089","GSM1807974","GSM1807990","GSM2055788","GSM2142335","GSM1807979","GSM1695909","GSM1554468","GSM1695910","GSM1807975","GSM1807977","GSM1807985","GSM1417183","GSM1807992","GSM1505572","GSM1807991","GSM1807996","GSM1416804","GSM2055782","GSM1807988","GSM1807981","GSM1807993","GSM1505598","GSM1807984","GSM1807973","GSM1807987","GSM1417184","GSM1807983","GSM1807982","GSM1807976",
"GSM1807980","GSM1807995","GSM1807994","GSM1807978","GSM2326088","GSM1807989","GSM1622386","GSM2199518","GSM1548454","GSM1707675","GSM2113310","GSM1548461","GSM1548457","GSM1548466","GSM1548460","GSM1548459","GSM1662561","GSM2199516","GSM1548455","GSM1662562","GSM2199514","GSM2199513","GSM1662703","GSM2199520","GSM2199519","GSM1662560","GSM2127797","GSM2199511","GSM1548463","GSM2113305",
"GSM2199512","GSM1662559","GSM1548465","GSM2199509","GSM2127795","GSM1960355","GSM1707674","GSM1548456","GSM1548453","GSM2127796","GSM1662558","GSM2199517","GSM1548462","GSM1662702","GSM1960357","GSM1548458","GSM1548464","GSM2199510","GSM1960356","GSM2199515","GSM1417182","GSM1505571","GSM1554467","GSM1807986","GSM2157100","GSM2157098","GSM2157102","GSM2388500","GSM2388502","GSM2388505",
"GSM2388508","GSM2388510","GSM2388501","GSM2388504","GSM2388506","GSM2388507","GSM2388503","GSM2388509","GSM2072372","GSM2072373","GSM2343657","GSM2343447","GSM2072387","GSM2072386","GSM2343656","GSM2343776","GSM2071276","GSM1067795","GSM1538102","GSM1538070","GSM1538073","GSM1538074","GSM1538079","GSM1538078","GSM1538098","GSM1538081","GSM1538080","GSM1538088","GSM1538093","GSM1538091",
"GSM1538069","GSM1538092","GSM1538101","GSM1538090","GSM1538094","GSM1538095","GSM1538089","GSM1901151","GSM1901147","GSM1901153","GSM1901149","GSM1974235","GSM1974236","GSM2042243","GSM2042242","GSM2170953","GSM2170940","GSM2170954","GSM2170939","GSM2170950","GSM2170957","GSM2170938","GSM2170958","GSM2170944","GSM2171016","GSM2170948","GSM2170949","GSM2170931","GSM2170946","GSM2170963",
"GSM2170966","GSM2170951","GSM2170952","GSM2170947","GSM2170965","GSM2170932","GSM2170935","GSM2170945","GSM2170962","GSM2171025","GSM2170928","GSM2171014","GSM2170961","GSM2170982","GSM2170964","GSM2170930","GSM2170936","GSM2170937","GSM2170923","GSM2170934","GSM2170879","GSM2170870","GSM2171022","GSM2170960","GSM2170924","GSM2170927","GSM2170874","GSM2170877","GSM2170926","GSM2170929",
"GSM2170987","GSM2171026","GSM2170925","GSM2170933","GSM2170955","GSM2170986","GSM2170878","GSM2170872","GSM2170988","GSM2170984","GSM2171015","GSM2170989","GSM2170983","GSM2170959","GSM2170993","GSM2171023","GSM2170977","GSM2170981","GSM2170922","GSM2171020","GSM2170883","GSM2170996","GSM2170868","GSM2170979","GSM2171024","GSM2171018","GSM2170903","GSM2170885","GSM2170890","GSM2170859",
"GSM2170876","GSM2171021","GSM2170969","GSM2170985","GSM2170956","GSM2170913","GSM2170999","GSM2170881","GSM2170898","GSM2170884","GSM2170991","GSM2170861","GSM2170860","GSM2170974","GSM2171004","GSM2170896","GSM2170899","GSM2170887","GSM2170889","GSM2170882","GSM2170875","GSM2170862","GSM2170869","GSM2171019","GSM2171017","GSM2171002","GSM2170915","GSM2170942","GSM2170888","GSM2170886",
"GSM2170873","GSM2170867","GSM2170855","GSM2170943","GSM2171009","GSM2170909","GSM2170976","GSM2170980","GSM2170894","GSM2170880","GSM2170865","GSM2170854","GSM2170895","GSM2170907","GSM2170900","GSM2170902","GSM2170905","GSM2170970","GSM2170864","GSM2170871","GSM2170911","GSM2170892","GSM2170972","GSM2170916","GSM2170904","GSM2170863","GSM2170891","GSM2170914","GSM2170893","GSM2171008",
"GSM2170856","GSM2170858","GSM2170897","GSM2170978","GSM2170866","GSM2170994","GSM2171005","GSM2170906","GSM2170918","GSM2171003","GSM2170995","GSM2170973","GSM2170857","GSM2170853","GSM2170908","GSM2170971","GSM2171013","GSM2170912","GSM2170975","GSM2170967","GSM2171001","GSM2170997","GSM2170968","GSM2170901","GSM2171000","GSM2170998","GSM2170910","GSM2171011","GSM2170941","GSM2171007",
"GSM2170919","GSM2170990","GSM2170920","GSM2170921","GSM2171010","GSM2171006","GSM2171012","GSM2170992","GSM2170917","GSM2232217","GSM2232200","GSM2232194","GSM2232199","GSM2232186","GSM2232220","GSM2232185","GSM2232204","GSM2232213","GSM2232182","GSM2232197","GSM2232219","GSM2232189","GSM2232198","GSM2232218","GSM2232223","GSM2232181","GSM2232209","GSM2232191","GSM2232212","GSM2232211",
"GSM2232179","GSM2232224","GSM2232214","GSM2232208","GSM2232216","GSM2232184","GSM2232221","GSM2232193","GSM2232210","GSM2232205","GSM2232195","GSM2232196","GSM2232203","GSM2232207","GSM2232180","GSM2232201","GSM2232187","GSM2232192","GSM2232188","GSM2232190","GSM2232202","GSM2232206","GSM2232183","GSM2232215","GSM2232222","GSM2262394","GSM2262403","GSM2262402","GSM2262397","GSM2262398",
"GSM2262406","GSM2262396","GSM2262401","GSM2262399","GSM2262400","GSM2262395","GSM2262405","GSM2262404","GSM2388499","GSM2388498","GSM2559748","GSM2559742","GSM2559741","GSM2559745","GSM2559739","GSM2559753","GSM2559752","GSM2559743","GSM2559749","GSM2559746","GSM2559750","GSM2559747","GSM2559751","GSM2559740","GSM2559755","GSM2559744","GSM2559738","GSM2559754","GSM2416873","GSM2416874",
"GSM2416875","GSM2416876","GSM2416877","GSM2416878","GSM2416879","GSM2416880","GSM2416881","GSM2416882","GSM2416883","GSM2416884","GSM2453421","GSM2610521","GSM2610522","GSM2610523","GSM2610524","GSM2753380","GSM2753381","GSM2753382","GSM2809205","GSM2809206","GSM2809207","GSM2809208","GSM2809209","GSM2809210","GSM2809211","GSM2809212","GSM2809213","GSM2809214","GSM2817112","GSM2453420",
"GSM2547178","GSM2547179","GSM2547180","GSM2547181","GSM2547182","GSM2547183","GSM2547184","GSM2547185","GSM2547186","GSM2547187","GSM2547188","GSM2547189","GSM2547190","GSM2547191","GSM2547192","GSM2547193","GSM2547194","GSM2547196","GSM2547197","GSM2547198","GSM2547199","GSM2547200","GSM2547201","GSM2547202","GSM2547203","GSM2547204","GSM2547205","GSM2547206","GSM2547207","GSM2547208",
"GSM2547209","GSM2547210","GSM2547211","GSM2547212","GSM2547213","GSM2547214","GSM2547215","GSM2547216","GSM2547217","GSM2547218","GSM2547219","GSM2547220","GSM2547221","GSM2547222","GSM2547223","GSM2547225","GSM2547226","GSM2547227","GSM2547228","GSM2547229","GSM2547230","GSM2547231","GSM2547232","GSM2547233","GSM2547234","GSM2547235","GSM2547236","GSM2547237","GSM2547238","GSM2547239",
"GSM2547240","GSM2547242","GSM2547243","GSM2547244","GSM2547245","GSM2547246","GSM2547247","GSM2547248","GSM2547249","GSM2547250","GSM2547251","GSM2547252","GSM2547253","GSM2547254","GSM2547255","GSM2547256","GSM2547257","GSM2547258","GSM2547259","GSM2547260","GSM2547261","GSM2547262","GSM2547263","GSM2547264","GSM2547265","GSM2547266","GSM2547267","GSM2547268","GSM2547269","GSM2547270",
"GSM2547272","GSM2547273","GSM2547274","GSM2547275","GSM2547276","GSM2547277","GSM2547278","GSM2547279","GSM2547280","GSM2547281","GSM2547282","GSM2547283","GSM2547284","GSM2547285","GSM2547286","GSM2547287","GSM2547288","GSM2547289","GSM2547290","GSM2547291","GSM2547292","GSM2547293","GSM2547294","GSM2547295","GSM2547296","GSM2547297","GSM2547298","GSM2547299","GSM2547300","GSM2547301",
"GSM2547302","GSM2547303","GSM2547304","GSM2547305","GSM2547306","GSM2547307","GSM2547308","GSM2547309","GSM2547310","GSM2547311","GSM2547312","GSM2547313","GSM2547314","GSM2547315","GSM2547316","GSM2547317","GSM2547318","GSM2547319","GSM2547320","GSM2547321","GSM2547322","GSM2547323","GSM2547324","GSM2547325","GSM2547326","GSM2547327","GSM2547328","GSM2547329","GSM2547330","GSM2547331",
"GSM2547332","GSM2547333","GSM2547334","GSM2547335","GSM2547336","GSM2547337","GSM2547338","GSM2547339","GSM2547340","GSM2547341","GSM2547342","GSM2547343","GSM2547344","GSM2547345","GSM2547346","GSM2547347","GSM2547348","GSM2547349","GSM2547350","GSM2547351","GSM2547352","GSM2547353","GSM2547354","GSM2547355","GSM2547356","GSM2547357","GSM2547358","GSM2547359","GSM2547360","GSM2547361",
"GSM2547362","GSM2547363","GSM2547364","GSM2547365","GSM2547366","GSM2547367","GSM2547368","GSM2547369","GSM2547370","GSM2547371","GSM2547372","GSM2547373","GSM2547374","GSM2547375","GSM2547376","GSM2547377","GSM2547378","GSM2547379","GSM2547380","GSM2547381","GSM2547382","GSM2547383","GSM2547384","GSM2547385","GSM2547386","GSM2547387","GSM2547388","GSM2547389","GSM2547390","GSM2547391",
"GSM2547392","GSM2547393","GSM2547394","GSM2547395","GSM2547396","GSM2547397","GSM2547398","GSM2547399","GSM2547400","GSM2547401","GSM2547402","GSM2547403","GSM2547404","GSM2547405","GSM2547406","GSM2547407","GSM2547408","GSM2547409","GSM2547410","GSM2547411","GSM2547412","GSM2547413","GSM2547414","GSM2547415","GSM2547416","GSM2547417","GSM2547418","GSM2547419","GSM2547420","GSM2547421",
"GSM2547422","GSM2547423","GSM2547424","GSM2547425","GSM2547426","GSM2547427","GSM2547428","GSM2547429","GSM2547430","GSM2547431","GSM2547432","GSM2547433","GSM2547434","GSM2547435","GSM2547436","GSM2547437","GSM2547438","GSM2547439","GSM2547440","GSM2547441","GSM2547442","GSM2547443","GSM2547444","GSM2547445","GSM2547446","GSM2547447","GSM2547448","GSM2547449","GSM2547450","GSM2547451",
"GSM2547452","GSM2547453","GSM2547454","GSM2547455","GSM2547456","GSM2547457","GSM2547458","GSM2547459","GSM2547460","GSM2547461","GSM2547462","GSM2547463","GSM2547464","GSM2547465","GSM2547466","GSM2547467","GSM2547468","GSM2547469","GSM2547470","GSM2547471","GSM2547472","GSM2547473","GSM2547474","GSM2547475","GSM2547476","GSM2547477","GSM2547478","GSM2547479","GSM2547480","GSM2563043",
"GSM2563044","GSM2677879","GSM2677880","GSM2677881","GSM2677882","GSM2677883","GSM2677884","GSM2701558","GSM2701559","GSM2701560","GSM2701561","GSM2701562","GSM2730123","GSM2730124","GSM2730125","GSM2730126","GSM2730127","GSM2730128","GSM2730129","GSM2730130","GSM2794842","GSM2794843","GSM2794844","GSM2794845","GSM2794846","GSM2794847","GSM2794848","GSM2794849","GSM2794850","GSM2794851",
"GSM2794852","GSM2794853","GSM2794854","GSM2794855","GSM2794856","GSM2794857","GSM2794858","GSM2794859","GSM2794860","GSM2794861","GSM2807391","GSM2807393","GSM2807395","GSM2807396","GSM2807397","GSM2807401","GSM2807403","GSM2807405","GSM2807407","GSM2807409","GSM2807411","GSM2807413","GSM2807415","GSM2807417","GSM2807419","GSM2807421","GSM2807422","GSM2807423","GSM2807427","GSM2807429",
"GSM2807431","GSM2807433","GSM2807435","GSM2807437","GSM2807438","GSM2807441","GSM2807443","GSM2807445","GSM2807447","GSM2807449","GSM2459128","GSM2459129","GSM2459130","GSM2459131","GSM2459132","GSM2459133","GSM2459134","GSM2459135","GSM2459136","GSM2459137","GSM2459138","GSM2459139","GSM2875111","GSM2875112","GSM2875113","GSM2875114","GSM2875115","GSM2875116","GSM2875117","GSM2875118",
"GSM2875119","GSM2875120","GSM2875121","GSM2875122","GSM2875123","GSM2875124","GSM2224921","GSM2224922","GSM2224923","GSM2224924","GSM2224925","GSM2224926","GSM2224927","GSM2224928","GSM2224929","GSM2224930","GSM2224931","GSM2224932","GSM2330093","GSM2330094","GSM2330098","GSM2330100","GSM2343909","GSM2343960","GSM2836481","GSM2480437","GSM2480439","GSM2480441","GSM2480443","GSM2480445",
"GSM2480447","GSM2480449","GSM2480451","GSM2480453","GSM2480455","GSM2480457","GSM2480459","GSM2480461","GSM2480463","GSM2480465","GSM2480467","GSM2480469","GSM2480471","GSM2480473","GSM2480475","GSM2480477","GSM2650363","GSM2650364","GSM2650365","GSM2650366","GSM2650367","GSM2650368","GSM2586505","GSM2586506","GSM2586507","GSM2586509","GSM2644719","GSM2644720","GSM2715277","GSM2715278",
"GSM1849373","GSM2653819","GSM2653820","GSM2653821","GSM2653822","GSM2653823","GSM2653825","GSM2653826","GSM2653827","GSM2653828","GSM2653829","GSM2653830","GSM2653831","GSM2653832","GSM2653833","GSM2653834","GSM2653835","GSM2653836","GSM2653837","GSM2653838","GSM2653840","GSM2653841","GSM2653842","GSM2653843","GSM2653844","GSM2653845","GSM2653846","GSM2653847","GSM2653849","GSM2653850",
"GSM2653851","GSM2653853","")

    # Retrieve information from compressed data
    samples = h5read(destination_file, "meta/Sample_geo_accession")
    tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
    genes = h5read(destination_file, "meta/genes")

    # Identify columns to be extracted
    sample_locations = which(samples %in% samp)

    # extract gene expression from compressed data
    expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
    H5close()
    rownames(expression) = genes
    colnames(expression) = samples[sample_locations]

    # Print file
    write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
    print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))
}
