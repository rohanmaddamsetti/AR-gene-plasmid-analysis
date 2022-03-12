#!/usr/bin/env python

'''
annotate-ecological-category.py by Vincent Huang and Rohan Maddamsetti.

This script annotates the ecology for each analyzed genome, based on the metadata in the
'host' and 'isolation_source' tags in its genome file from Genbank.

Some definitions for annotation: soil is everything from dirt, whether an agricultural field or not.
agriculture, on the other hand, refers strictly to crop plants-- soil is not included.

In addition, "rhizosphere" always maps to the Soil annotation.

"Animal" also includes fungi and protists-- basically all non-plant eukaryotes.

Usage: python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv
'''

########################################################################################################################
## Keywords for Annotating Ecological Categories.
########################################################################################################################

animal_hosts = [
    "paulinella",
    "gut of protaetia brevitarsis seulensis (larva)", "entylia carinata", "ectobius sp.",
    "oncorhynchus mykiss (Rainbow trout)", "drosophila mauritania", "meriones meridianus",
    "ctenopharyngodon idella (grass carp)", "siniperca chuatsi", "oncorhynchus tshawytscha",
    "tabanus nigrovittatus (horsefly)", "reticulitermes speratus", "curculio caryae",
    "rabbit tick", "rhopalodia gibberula", "ornithodoros parkeri",
    "auritibicen bihamatus", "prairie dog", "protaetia brevitarsis seulensis larva",
    "diestrammena coreana (camel cricket)", "sycon capricorn", "euscepes postfasciatus",
    "mediastinia sp.", "cosmozosteria sp.", "amblyomma americanum", "platyzosteria sp.",
    "cellana toreuma", "ring-tailed lemur", "testudinidae", "photinus",
    "haemaphysalis concinna", "apis cerana", "heterodera glycines", "pyxicephalus edulis",
    "trialeurodes vaporariorum", "ciconia boyciana (oriental stork)", "tettigades chilensis",
    "cryptotympana atrata", "bemisia tabaci china 1", "kentomonas sorsogonicus",
    "animal", "lizard", "wildebeest", "urostylis westwoodii", "steinernema monticolum",
    "phoxinus keumkang", "sturgeon", "heterorhabditis bacteriophora", "gerbil",
    "ixodes ovatus", "otter", "odobenus sp.", "pachyrhynchus infernalis",
    "epinephelus aeneus", "sander vitreus", "sheep (ovis)", "anser albifrons",
    "pantholops hodgsonii", "cats", "capra aegagrus hircus", "cimex lectularius jesc",
    "spermophilus dauricus", "diceroprocta semicincta (cicada)", "frog",
    "donacia thalassina", "ixodes persulcatus", "tegillarca granosa", "wild bird",
    "edessa sp.","nematode", "seriola lalandi", "marmota himalayana", "egret",
    "tragelaphus strepsiceros", "apis mellifera", "botryllus sp.", "oncorhynchus kisutch",
    "meleagris gallopavo", "blue sheep", "shellfish", "insect", "chelymorpha alternans",
    "neophocaena asiaeorientalis", "penaeus japonicus", "alligator", "macroplea mutica",
    "tibetan antelope", "penaeus vannamei (whiteleg shrimp)", "larus sp.", "protozoa",
    "donacia provostii","donacia marginata", "ictalurus punctatus", "amblyomma cajennense",
    "argas persicus","mink", "cistudinella sp.", "ischnocodia annulus", "bigeye tuna",
    "drosophila melanogaster oregon-r-modencode", "sablefish", "oncorhynchus mykiss",
    "atlantic salmon", "shrimp", "amblyomma variegatum (cattle tick)", "nile tilapia",
    "tilapia fish", "migratory bird", "south china tiger", "crocodile lizard",
    "arvelius albopunctatus", "crassostrea gigas","phoca largha", "misgurnus anguillicaudatus",
    "xenopus laevis", "trigona sp.", "buteo jamaicensis", "aphis fabae", "eclectus roratus",
    "equus kiang", "phractocephalus hemioliopterus", "black-collared starling", "marmot",
    "cryptocercus clevelandi", "wild chukar", "alvinocaris longirostris",
    "penaeus vannamei","donacia versicolorea","nipponaphis monzeni","salmo salar",
    "amazona sp.", "periplaneta americana", "mizuhopecten yessoensis", "haliotis discus hannai",
    "plateau pika", "silkworm", "gentoo penguin", "discomorpha sp.", "tibetan gazelle",
    "chroicocephalus novaehollandiae (australian silver gull chick)", "trachinotus ovatus",
    "crow", "cassida sp.", "pelodiscus sinensis", "hirudo verbana", "procambus clarkii",
    "galaxea fascicularis (stony coral)","bivalve mollusk","ondara zibethicus",
    "aphis glycines", "mus musculus", "sebastes schlegeli", "mus musculus subsp. domesticus",
    "donacia semicuprea","tuberolachnus salignus", "hippopotamus amphibius", "jellyfish",
    "paralichthys olivaceus", "bos mutus", "pheasant", "charidotella sexpunctata",
    "scylla serrata", "euphausia superba", "canis", "adelges kitamiensis",
    "neoaliturus tenellus", "crucian carp", "ictalurus furcatus", "mercenaria mercenaria",
    "physopelta gutta", "agroiconota sp.", "pentalonia nigronervosa", "bison bison",
    "apterostigma dentigerum","urocitellus undulatus", "mussels", "artemisaphis artemisicola",
    "aphis nerii","eel", "ixodes ricinus (tick)", "acromis sparsa", "carposina sasakii",
    "macroplea appendiculata","donacia clavipes","connochaetes taurinus","perinereis linea",
    "trichonympha agilis", "wild boar", "pseudotrichonympha sp.", "flounder",
    "drosophila melanogaster","dicentrarchus labrax (sea bass)", "rainbow trout",
    "donacia dentata", "oyster larvae", "myocastor coypus", "megalobrama amblycephala",
    "honey bee", "channa argus (snakehead fish)", "crocodylus siamensis", "forficula auricularia",
    "ixodes scapularis","gut of wasp", "american cockroach", "trypoxylus dichotomus",
    "parrot", "paper wasp", "marmota", "scophthalmus maximus", "locusta migratoria",
    "oreochromis niloticus", "ornithodoros hermsi", "oecophylla smaragdina", "euschistus servus",
    "ixodes ricinus", "cassida viridis", "cassida vibex", "cynoglossus semilaevis",
    "phakellia ventilabrum", "phoca vitulina", "megacopta punctatissima", "invertebrates",
    "acyrthosiphon pisum", "donacia vulgaris","mouse", "mya arenaria oonogai makiyama",
    "nasonia vitripennis","cryptocercus punctulatus", "porcellio scaber", "oreochromis sp.",
    "schizaphis graminum biotype i", "oreochromis", "cage-cultured red drum", "coral",
    "ornithodoros sonrai", "fish", "amblyomma longirostre", "ceratitis capitata",
    "gorilla", "hydra vulgaris aep", "amblyomma neumanni","plateumaris rustica","oryctes gigas",
    "crab","honeybees", "aratinga solstitialis", "correlophus ciliatus", "struthio camelus",
    "syngnathus typhle","nezara viridula (cotton pathogen vector; southern green stink bug)",
    "deinagkistrodon acutus", "aphis helianthi", "anguilla anguilla", "halyomorpha halys",
    "antho dichotoma", "grass carp", "bighead carp", "trichechus manatus",
    "pthirus gorillae","marmota baibacina","ruddy shelduck", "marmota sibirica",
    "protaetia brevitarsis seulensis", "cancer pagurus", "lepus brachyurus", "lepus europaeus",
    "bactrocera oleae","rattus", "hyadaphis tataricae", "myzus persicae", "alces alces",
    "flatfish","toothfish", "bankia setacea", "litopenaeus vannamei","naegleria",
    "acyrthosiphon kondoi", "sebastes schlegelii", "penaeus vannamei (shrimp)", "ovis",
    "macrosiphum gaurae", "turbot", "turbo", "adelges lariciatus", "adelges cooleyi",
    "pseudotrichonympha grassii (protist) in the gut of the termite coptotermes formosanus",
    "cerambycidae sp.", "pediculus humanus corporis", "ixodes pacificus","donacia proxima",
    "kangaroo", "cryptopsaras couesii","anas platyrhynchos", "penaeus setiferus",
    "apis mellifera mellifera", "culicoides impunctatus", "larimichthys crocea",
    "camponotus chromaiodes", "caenorhabditis elegans","murgantia histrionica",
    "parakeet", "geronticus eremita","melanaphis sacchari","dendrolimus ibiricus",
    "drosophila melanogaster oregon", "allacta bimaculata", "peafowl", "melopsittacus undulatus",
    "homalodisca coagulata (glassy-winged sharpshooter)", "pachypsylla venusta",
    "myzus persicae (green peach aphid)", "macrotermes barneyi", "glossina morsitans morsitans",
    "nauphoeta cinerea", "delphinapterus leucas (beluga whale)", "sea urchin",
    "oryctolagus cuniculus","corvus brachyrhynchos","drosophila melanogaster oregon-r modencode",
    "apterostigma", "helicoverpa armigera", "seabass", "rabbit", "plateumaris braccata","rat",
    "odocoileus virginianus","peromyscus leucopus", "cinara cedri", "brachycaudus cardui",
    "plateumaris consimilis", "acanthamoeba polyphaga hn-3", "perameles bougainville",
    "hyalomma aegyptium","zebra", "apsterostigma", "baizongia pistaciae", "pediculus schaeffi",
    "otospermophilus beecheyi (california ground squirrel)","apodemus agrarius",
    "paralichthys olivaceus (flounder)","dermacentor andersoni","donacia tomentosa","lacertilia",
    "papio papio", "chroicocephalus novaehollandiae", "ochotona curzoniae (plateau pika)",
    "anomalops katoptron", "gull", "japanese eel", "salvelinus fontinalis","blattella germanica",
    "pteropus livingstonii", "pteropus poliocephalus", "bat", "skate", "red kangaroo",
    "moschus berezovskii", "macaca silenus", "panulirus ornatus","pavona duerdeni",
    "odontobutis platycephala", "lutjanus guttatus (rose snapper)", "kentomonas sorsogonicus",
    "ostrea edulis (flat oyster)", "beluga whale", "shinkaia crosnieri",
    "pelteobagrus fulvidraco", "aquila chrysaetos", "graptopsaltria bimaculata",
    "vultur gryphus","penaeus monodon","tanakia koreensis", "dicentrarchus labrax",
    "zophobas atratus", "pigeon", "turtle", "penaeus (litopenaeus) vannamei (whiteleg shrimp)",
    "donacia cinerea", "fulmars", "brevicoryne brassicae", "acanthamoeba", "blaberus giganteus",
    "arctic char", "microlophium carnosum", "photinus pyralis", "ciconia boyciana",
    "macaca mulatta", "thelaxes californica","bactrocera dorsalis","apostichopus japonicus",
    "ovis aries", "leptinotarsa decemlineata", "phyllophaga sp.","bombyx mori",
    "pthirus pubis","manis javanica (pangolin)", "sipalinus gigas", "geodia barretti",
    "macrosiphoniella sanborni", "muscaphis stroyani","aphis urticata",
    "crassostrea gigas (pacific oyster)", "crassostrea virginica", "ellychnia corrusca",
    "glossina brevipalpis","marine sponge lissodendoryx isodictyalis in the bahamas",
    "catfish", "morone chrysops x morone saxatilis", "pediculus humanus", "rattus norvegicus",
    "bothriocroton concolor", "hyperomyzus lactucae","trimyema compressum", "macaca fascicularis",
    "neohaemonia nigricornis", "haemaphysalis juxtakochi", "columba livia", "diaphorina citri",
    "prawn","deer", "nebria ingens riversi", "donacia fulgens", "sitobion avenae",
    "pseudorca crassidens","acyrthosiphon lactucae","aphis nasturtii","giant panda",
    "cryptocercus kyebangensis", "donacia sparganii", "oyster", "caprine",
    "japanese rhinoceros beetle larva","hydrophilus acuminatus", "discus", "hawk",
    "acinonyx jubatus","chaeturichthys stigmatias", "donacia cincticornis", "anas strepera",
    "pan troglodytes", "paguma larvata", "plateumaris sericea", "diplonemea",
    "aphis craccivora (cowpea aphid)", "calanoid copepod", "atlantic white-sided dolphin",
    "ctenocephalides felis", "parachirida sp.", "nezara viridula", "euterpnosia chibensis",
    "pyropia yezoensis conchocelis", "pan paniscus", "lipaphis pseudobrassicae",
    "donacia piscatrix", "diaphorina cf. continua", "lion-tailed macaques", "penguin",
    "aphis craccivora","donacia crassipes", "donacia bicoloricornis", "donacia simplex",
    "draeculacephala minerva", "spodoptera frugiperda", "danaus plexippus", "danio rerio",
    "anodonta arcaeformis", "mastotermes darwiniensis", "toxic alexandrium minutum",
    "coreoleuciscus splendidus", "avian", "mus musculus c57bl/6j", "tick",
    "platycercus elegans", "solea senegalensis","sanzinia madagascariensis volontany",
    "plateumaris pusilla","gilthead seabream", "ixodes pacificus (western blackleg tick)",
    "coregonus clupeaformis (lake whitefish)", "wax moth", "galleria mellonella",
    "cyclopterus lumpus","cockatiel", "melanocetus johnsonii", "acheta domesticus",
    "blatta orientalis", "seriola dumerili", "white stork", "canis latrans",
    "drosophila neotestacea","shrimps", "bird of prey", "ardea cinerea", "korean rockfish",
    "c57bl/6ntac mice", "phascolarctos cinereus", "neopsylla setosa", "haliclona simulans",
    "rhopalosiphum padi", "silurus asotus", "uroleucon ambrosiae","mustela putorius furo",
    "meimuna kuroiwae", "tettigades ulnaria", "cystophora cristata", "porites pukoensis",
    "ringed seal", "holothuroidea", "tilapia", "macropus giganteus",
    "exterior surface of the shell of an abalone sold in a fish market in tokyo; japan",
    "ailuropoda melanoleuca","spermophilus sp.", "squirrel", "dendrolimus sibiricus", "cormorant",
    "melanozosteria sp.", "protagonista lugubris", "therea regularis"]

animal_isolation_sources = [
    "isolated from haemaphysalis hystricis",
    "feces of an oriental stork (ciconia boyciana)", "reptile",
    "crayfish", "iguana iguana", "crustacean shell", "malayan pangolin",
    "animal", "wild rat", "coho salmon", "fruit fly","diaphorina citri", "south china tiger",
    "iguana", "bat", "wild yak feces", "white-crowned sparrow", "african elephants",
    "efb-infected honey bee colony", "turban shell", "mussels", "rainbow trout",
    "foulbrood of honeybees", "turbot", "red fox", "panda", "neopsylla setosa",
    "marine sponge", "chamaeleonidae","ixodes spinipalpis", "barnacle at wood pile-on",
    "atlantic salmon", "manis javanica", "a feces sample of migratory birds origin",
    "tilapia", "wild pig; fecal", "shellfish","diarrheal snake diarrheal snake in hunan",
    "diseased shrimp", "wasp honeycombs", "eclectus roratus feces", "frog",
    "insect larvae", "tissue; animal","yellowtail", "snake", "hydrophilus acuminatus",
    "crassostrea gigas","pigeon", "gill tissue of bathymodiolus japonicus",
    "meleagris gallopavo","sea bass","mouse gut","avian","sockeye salmon",
    "diseased labeo rohita fish","mealworm", "feces; gull (larus spp.)", "nematode community",
    "musca domestica","atlantic cod","siniperca scherzeri", "red-breasted parakeet",
    "midgut crypts of stink bug togo hemipterus","konosirus punctatus","heterodon nasicus",
    "crushed cell of callyspongia sp.", "tuna", "housefly", "bison", "coral","catfish",
    "mactra veneriformis", "black-collared starling", "white-lipped deer", "elk droppings",
    "shrimp", "black-headed gull", "psittacus erithacus feces", "salmon", "silkworm excrement",
    "separated from the corpses of silkworms that had died due to bb natural infection in daiyue district; taian city; shandong province; china.", "rabbit feces", "jellyfish",
    "rabbit", "animal hide","tuna scrape; yellowfin", "myocastor coypus",
    "cybister lewisianus", "intestinal contents of plateau pika", "vole",
    "tuna scrape; frozen yellowfin nakaoachi", "brook charr", "penaeus japonicus",
    "aphis craccivora (cowpea aphid) on robinia pseudoacacia (locust)", "insects",
    "swiss alpine ibex feces", "psittacus erithacus", "penaeus vannamei (whiteleg shrimp)",
    "bearded dragon", "wild bird", "rattus rattus", "hepatopancreas", "pagrus major",
    "neophocaena phocaenoides", "loach","unidentified actinians", "fish",
    "waxworms gut", "flea", "intestinal contents of termite nasutitermes nigriceps",
    "muskrat during outbreak", "invertebrate gut", "cockless", "herring gull; cloacal swab",
    "gromphadorhina portentosa cockroaches"]

plant_hosts = [
    "suaeda salsa","cladonia borealis","medicago truncatula", "ulva pertusa (algae)",
    "allamanda cathartica","astragalus pelecinus","lebeckia ambigua","lotus corniculatus",
    "trifolium repens", "curcuma aromatica", "curcuma wenyujin y.h. chen et c. ling",
    "musa balbisiana cultivar kepok", "lemna trisulca", "pelargonium capitatum",
    "datisca glomerata","beta vulgaris","robinia pseudoacacia","ficus religiosa l.",
    "sporobolus anglicus", "casuarina equisetifolia", "mulberry", "ficus benjamina",
    "leontopodium alpinum","tanacetum vulgare","catharanthus roseus","stereocaulon sp.",
    "oxytropis pumilio", "malus sylvestris", "angelica sinensis dlies",
    "agave americana l.","stachytarpheta glabra","elymus tsukushiensis","eucalyptus",
    "androsace koso-poljanskii", "bromus inermis", "peperomia dindygulensis miq.",
    "trifolium pratense","lemna minor","mimosa affinis","salsola stocksii", "mimosa",
    "biserrula pelecinus l.","kaempferia sp. rhizome","masson pine", "populus x jackii",
    "nerium oleander","biserrula pelecinus","lichen","miscanthus giganteus","mimosa scabrella",
    "rosa sp.", "malus sieversii", "onobrychis viciifolia", "carex sp. (sedge blades)",
    "flower","mimosa flocculosa","trifolium uniflorum","eucalypti of eucalyptus","wild cotton",
    "urochloa reptans", "cotinus coggygria", "kandelia candel (mangrove)",
    "trifolium spumosum l. (annual mediterranean clovers)","microcystis culture",
    "medicago orbicularis","arabidopsis", "vachellia farnesiana", "arabidopsis thaliana",
    "physcomitrella patens","cryptomeria fortunei","festuca arundinacea","oxytropis triphylla",
    "usnea","rubus sp.", "melilotus officinalis", "himalayan blackberry",
    "achillea ptarmica l.","muscari","lespedeza cuneata","nerium oleander",
    "tanacetum vulgare","eucommia ulmoides", "clematis", "quercus rubra",
    "leiosporoceros dussii","oxytropis triphylla","vavilovia formosa", "oxytropis kamtschatica",
    "zantedeschia aethiopica","rosa", "prunus cerasifera", "medicago",
    "leontopodium nivale","nusuttodinium aeruginosum","prosopis cineraria",
    "oryza glumipatula","acacia farnesiana", "catalpa", "populus", "jatropha curcas",
    "populus alba x (p. davidiana + p. simonii) x p. tomentosa",
    "salvia splendens","panax ginseng", "vicia alpestris","phaeoceros",
    "commiphora wightii", "alhagi sparsifolia shap.", "acacia farnesiana",
    "lobaria pulmonaria thallus", "lotus", "erigeron annuus l. pers",
    "hyacinthus orientalis","sporobolus anglicus","malus prunifolia (crab apple)",
    "tephrosia apollinea", "vitis", "echinacea purpurea", "phyla canescens",
    "aeschynomene fluminensis","vavilovia formosa", "prunus sp.",
    "oenothera speciosa", "acaciella angustissima", "brassica rapa subsp. chinensis",
    "catharanthus roseus","dracaena sanderiana","parthenium argentatum gray (guayule shrubs)",
    "medicago arborea", "medicago lupulina l.", "dongxiang wild rice",
    "euonymus sp.","eucalyptus grandis","broussonetia papyrifera","lichen stereocaulon sp.",
    "lotus sp.","weed", "thlaspi arvense", "grass",
    "blasia pusilla","calliandra grandiflora","brittle root","digitaria eriantha","vavilovia formosa",
    "plant", "sesame", "salix sp.", "hibiscus", "hibiscus rosa-sinensis"]

plant_isolation_sources = [
    "flower of forsythia koreana", "algal phycosphere", "lotus corniculatus nodule",
    "malus sylvestris", "morning glory", "dongxiang wild rice", "pistachios",
    "flower of rhododendron schlippenbachii", "root of codonopsis pilosula", "root nodules",
    "chrysochromulina tobin phycosphere", "algae","artificially infected catharanthus",
    "leaf spots", "poplar tree gall", "gardenia fruit", "plant root",
    "peltigera membranacea thallus","angelica gigas nakai root surface",
    "galega orientalis root nodule", "seaweed", "xylem", "red seaweed",
    "plant","broussonetia papyrifera","lichen", "lowbush blueberries",
    "seaweed; enteromorpha linza (coastal marine sediments); aburatsubo inlet; australia",
    "flower","leaf of fig tree", "grass from the park in a vetenary clinic",
    "plant roots","galega officinalis root nodule","lichen thallus","ginseng","gladiolus",
    "sesbania spp. root nodule","plant leaf","origanum marjorana","malus sp.","algae",
    "flower of shittah tree", "himalayan blackberry",
    "inner tissues of halophyte limonium sinense (girard) kuntze","leaf",
    "root nodule of sesbania cannabina", "wood",
    "fallen leaves of virgin forest", "flower of chiness redbud",
    "flower of rhododendron sclippenbachii", "plant xylem"]

terrestrial_isolation_sources = [
    "solfataric thermal field close to moutnovsky volcano",
    "iron hydroxide deposits", "solar salt", "air", "elton hypersaline lake",
    "brazilian saline-alkaline lake", "desert", "saltpan", "hyperthermophilic compost",
    "coal bed","moderate hot spring","salt marsh","alkaline pool submerged anode electrode",
    "kulunda steppe hypersaline lake", "hypersaline environment",
    "chronically low temperature and dry polar region", "thermophilic environment",
    "solar saltern", "dust", "beach sand", "salar de atacama; atacama desert",
    "weathered rock sample", "antimony mine", "sand", "weathered tuff",
    "moist arsenopyrite (feass)-containing rock taken from a mine tunnel approximately 300 m below the ground in the granites gold mine", "microbial mat", "subsurface rock",
    "salt crystallizer of little rann of kutch", "stalactite biofilm",
    "inside the caves of drach", "non-purified solar salt","gold-copper mine",
    "xinjiang aibi salt lake", "carbonated precipitates", "soda-saline lake",
    "shar-burdin hypersaline soda lake", "marine solar saltern brine","soda lake magadi",
    "rock (ENVO:00001995)", "Geothermal Reservoir", "rocky sand", "solar saltern in gomso bay",
    "lava", "brine", "saltern", "baengnokdam summit crater area; mt. halla",
    "taibei marine solar saltern near lianyungang city","orthoquartzite cave surface","cave",
    "the bange salt-alkaline lake in tibet", "the surfaces of weathered potassic trachyte",
    "deep subseafloor coal bed", "deep subsurface anoxic brine",
    "phototrophic microbial mat in hot lake; a shallow mgso4 dominated salt lake",
    "geothermal isolate", "salt crust", "coal spoil heap",
    "flat; laminated microbial mat in a salt marsh", "solar saltern of 19% salinity",
    "rock (envo:00001995)", "saline saltern", "hypersaline lake", "yates shaft; surf",
    "thiodendron' bacterial sulfur mat from mineral sulfide spring"]

soil_isolation_sources = [
    "agricultural field", "tundra wetland", "permafrost", "compost", "mixed sand sample",
    "sphagnum peat from the bog obukhovskoe (acidic wetland)",
    "paddy field; sungai manik; malaysia", "dune grassland",
    "rice fields", "sphagnum peat", "soil around hot spring", "dirt",
    "sphagnum bog", "fertilizer", "ancient permafrost from mammoth","montane grasslands",
    "solar salt farm", "enriched culture of compost", "farmland", "permafrost; kolyma lowland",
    "siberian permafrost","particulate matter", "pine forest", "mushroom compost",
    "manure compost", "angelo meadow plot 1; 20cm depth; 2 days after second rain event (91mm)",
    "composted cattle manure", "agricultural waste material",
    "long-term organic manure fertilized",
    "Microscale soil grain", "ermafrost region of qilian mountains"]

agri_hosts = [
    "coffea", "grape", "citrus spp.", "prunus avium", "capsicum sp.","phaseolus vulgaris",
    "sugarcane", "lolium perenne", "zingiber officinale", "raspberry",
    "phaseolus", "punica granatum", "pyrus pyrifolia var. culta","alfalfa","cucumis sativus",
    "cucurbita maxima", "lentil", "lathyrus sativus", "coffea arabica",
    "solanum melongena", "vitis vinifera","vitis vinifera l. cv. seto giants", "grapevine",
    "sakura tree", "rice", "japanese radish", "black pepper plant", "radish",
    "chinese cabbage","kiwifruit", "pisum sativum l. (pea)", "common bean", "melon", "pepper",
    "glycine max", "glycine", "manihot esculenta", "capsicum annuum",
    "actinidia deliciosa", "sugarcane", "actinidia chinensis", "lettuce", "olea europaea",
    "grapefruit", "brassica juncea var. foliosa", "apple", "strawberry",
    "brassica oleracea var. botrytis", "camellia sinensis", "juglans regia", "coffee",
    "brassica rapa subsp. pekinensis", "humulus lupulus", "morus alba",
    "triticum aestivum", "areca catechu", "date palm", "zea mays l.", "ginger","plum", "cowpea",
    "potato", "wheat", "malus domestica", "brassica oleracea", "rice plant",
    "saccharum officinarum", "cotton", "soybean","mandarin orange", "fava bean",
    "vitis vinifera l. cv. aurora black","ziziphus mauritiana lam", "pyrus pyrifolia",
    "sesame seedling", "tobacco","peanut", "carrot", "apple tree", "banana",
    "soybean (glycine max (l.) merrill)", "ipomoea aquatica",
    "panicum miliaceum", "solanum lycopersicum", "plantain", "pear", "winter wheat",
    "pyrus communis 'williams'", "fragariae ananassa", "hordeum vulgare",
    "actinidia", "zea mays", "turfgrass", "prunus dulcis","glycine soja", "potato",
    "pepper plant", "sweet orange", "mangifera indica", "maize",
    "cicer arietinum","pisum sativum", "actinidia deliciosa 'hayward'", "arachis hypogaea",
    "forage rape", "pear", "vigna radiata", "corn", "triticum aestivum (aestivum group)",
    "prunus cerasus (sour cherry)", "prunus dulcis", "sorghum bicolor", "phaseolus sp.",
    "allium cepa", "brassica oleracea var. capitata",
    "vitis vinifera cv. 'izsaki sarfeher'","brassica rapa var. laciniifolia subvar. oblanceolata",
    "allium cepa l.", "blueberry", "pogostemon cablin",
    "gossypium hirsutum","solanum lycopersicoides", "vicia faba", "eggplant", "pear tree",
    "glycine max cv. jinju1", "lycopersicon esculentum", "camellia oleifera",
    "allium cepa (onion)", "brassica rapa ssp. pekinensis (chinese cabbage)",
    "musa sp.", "apium graveolens",
    "triticum aestivum l.", "coffee plant", "pineapple", "citrus sinensis", "gossypium sp.",
    "musa spp."]

agri_isolation_sources = [
    "styrian pumpkin anthrosphere", "citrus paradisi",
    "greenhouse grown plants with black leg infection", "marjoram", "affected cabbage tissue",
    "ucb-1 pistachio rootstock", "pineapple", "plants with moko disease","mangoes", "chives",
    "raw almonds", "chilli", "eggplant", "cucumber", "rice leaves", "campbell early grape",
    "isolated from cucumber", "cabbage","almond drupe", "mexican lime leaf", "apple twig",
    "phaseolus vulgaris root nodule", "cucumber 'kurazh f1'", "barley grains",
    "solanum tuberosum", "peaches", "cornstalks and leaves", "maize seeds", "quinoa roots",
    "tomato roots", "pistachios; raw", "infected wheat", "sesame leaf", "spinach",
    "streptomyces spp. which was isolated from pyeonchang", "wheat roots", "corn root",
    "healthy tomato plant","almond kernel (raw; variety carmel)", "the root of rice", "soybeans",
    "wheat germ", "raw peanuts", "rot potato tubers", "litchi pericarp", "rice fields",
    "sugarcane root", "originally isolated from olive trees", "rice seed", "tea rhizoplane",
    "isolated from the tobacco substrate", "endorhiza sugar beet", "orange tree", "rice shoot",
    "hydroponic pots with potatoes", "ragi", "canola roots","leaf from rome apple cultivar",
    "orange","maize stem", "strawberry leaf tissue", "rye silage", "fodder","stable grass silage",
    "beans","carrot", "potato", "bean blight", "wilting pepper (c. annuum) stems in sanya",
    "corn", "japanese pear","strawberry leaf tissue", "wheat root", "ogi (red sorghum)",
    "oryza sativa", "nodules from common bean", "pisum sativum root-nodule", "grass silage",
    "noduls from roots of medicago sativa", "cilantro","rice root", "rice leaf", "maize leaf",
    "naturally-infected soybean tissue", "soybean", "citrus orchard", "tomato",
    "naturally-infected soybean leaf tissue", "apple tree", "banana", "cabbage seeds",
    "white nectarines fire pearl variety", "soybean nodules", "pistachio", "field isolate"]

marine_hosts = ["argopecten purpuratus", "seawater", "pyropia tenera",
                "trichodesmium erythraeum ims101", "asterionellopsis glacialis strain a3",
                "rhodosorus marinus", "red alga", "pyropia", "galaxea fascicularis",
                "heterostera chilensis", "grateloupia sp.", "brown algae"]

marine_isolation_sources = [
    "edge of california current after nitrate enrichment and low light incubation",
    "black smoker wall; 3500 m depth", "surface of plastic from ocean",
    "salt lake", "hydrothermal vent area derived", "deep sea", "ocean water",
    "water from the baltic sea", "envo:00000227", "sea", "envo:00002149",
    "sea water", "deep-sea hydrothermal vent", "sea ice", "oceanic water",
    "cold seep", "envo:01000301", "tidal pool", "coral; primnoid",
    "ocean", "isolated from the deepest ocean", "surface of fucus serratus",
    "chrysochromulina tobin phycosphere", "gulf of finland", "cyanobacterial bloom",
    "glass slide place on reef flat in natural environment",
    "surface water of the southern north sea","antarctic iceberg",
    "hydrothermal vent", "acidic salty water", "marker113 at caldera of axial seamount",
    "shallow tropical waters; normally from coral reef substrate","open ocean water",
    "sulfidic waters (60 m) from the peruvian upwelling region","shellfish hatchery",
    "bay of bengal", "pacific plankton",
    "culture of roseovarius sp. tm1035", "coastal","runway 10 reef (10-12 m)",
    "german wadden sea", "surface of pikea pinnata",
    "microbial mat material from brackish estuary", "hydrothermal precipitates",
    "shallow sea; symbiosis with algae", "surface of macroalgae",
    "beach sand","estuarine water [envo:01000301]", "150 km offshore",
    "cyanobacterial aggregates","synechocystis sp. gt-l",
    "seawater", "shallow-sea hydrothermal system"]

freshwater_hosts = ["water (environment)", "water", "mesostigma viride"]

freshwater_isolation_sources = [
    "creek","pond","tianshan glacier", "environmental (pond)", "oligrotrophic pond",
    "glomma river","snow","microcystis culture","glacial stream",
    "envo:00002011", "stream biofilm",
    "euglena gracilis from city ponds", "glacier from lahual spiti valley","river",
    "pool at botanical garden",
    "ayakekumsalt lake", "envo:00002006"]

food_hosts = ["kimchi cabbage", "mung bean", "kombucha scoby", "seujeot", "pickle",
              "pickled cabbage", "kimchi", "milk", "carrageenan", "frozen raw chicken",
              "beef", "koumiss", "dairy cattle", "nostoc flagelliforme", "sichuan pickle"]

food_isolation_sources = [
    "saeng-gimol meju", "rubing cheese", "instant soup", "soybean paste",
    "futsai", "koumiss", "soybean paste (chonggugjang)", "mashed potatoes", "fermented food",
    "kombucha scoby", "kefir", "cheese product", "instant pork", "pork &cabbage dumplings",
    "isolated from commercially available caspian sea yogurt", "broth",
    "makgeolli (korean traditional alcoholic beverage)", "cream products",
    "light wheat beer", "sprouts", "3-day-old traditional semihard cheese","pork",
    "moto starter of sake", "minced pork", "shrimp paste", "tahini (ground sesame seeds)",
    "isolated from salted brown alga laminaria", "shima dofu (okinawan-style tofu)",
    "salted brown alga laminaria", "stinky xiancaigeng", "sake(hatsuzoe)",
    "blue cheese in wax", "tibicos symbiotic community", "pickled cabbage",
    "ready to eat mixed salad leaves (obtained from discount store)",
    "degraded sugar thick juice","fresh alfalfa sprouts", "scallop", "chinese sausages",
    "ground turkey", "gimbap", "product-eggs-raw-whites", "traditional dairy product",
    "frozen peas", "palm wine", "natural whey culture from gruyere cheese",
    "homemade koumiss", "pea soup; fp", "madre curry powder hot", "highland barley wine",
    "salad", "rotten eggs", "steamed conch", "rice wine rice syrup", "garam masala",
    "water kefir", "palm sap","peanut butter", "ground turkey", "sick cider", "onion",
    "austrian hard cheese rind", "litopenaeus vannamei purchased at the supermarket",
    "chinese pickle", "cheese", "raw mutton","nuruk; korean traditional beverage starter",
    "crab marinated in soy sauce", "pork product", "raw pistachios", "chinese sauerkraut",
    "honey", "ham", "bobby veal steak", "10' wieners", "natto", "yogurt",
    "gochujang", "cheonggukjang", "cheese starter culture", "corn steep liquor",
    "traditional greek kasseri cheese","pork & cabbage dumplings", "lettuce",
    "new zealand cheese", "wine", "pork casserole; fp", "the red brine of salted laminaria",
    "beer", "tibetan kefir", "carrot juice", "pickle", "butter starter",
    "product-eggs-raw-whole", "tomato pulp", "frozen whole tilapia", "fish ball",
    "raw sausage", "beer contaminant", "scoby from kombucha tea", "calf liver",
    "mustard pickles", "organic spinach", "smoked fish", "cooked ox blood", "baijiu",
    "doubanjiang", "retail turkey","fresh produce (lettuce; lollo bionda)", "blue berry",
    "apple cider","fried eel", "radish and carrot pickled with rice bran and salt",
    "mac and cheese", "shelled pistachios", "marinated fish product", "vinegar pei",
    "10 weeks old 45+ samso cheese", "salami", "chinese traditional sourdough",
    "mixed salads", "tibet kefir","sugar thick juice", "traditional dairy products",
    "ice cream", "water containing garlic and cabbage", "doenjang(soybean paste)",
    "green chile ingredient", "ground beef", "soy sauce mash", "yoghurt",
    "koumiss (fermented mare's milk)", "pilsner beer","blown cheese",
    "pickled green chili peppers", "fish part (slab) sold as silver carp; muscle",
    "wheat sourdough", "cherry","pork barbeque", "vinegar", "dried tofu",
    "sauerkraut", "diced lamb", "salt", "commercial salt", "patagonian red wine",
    "korean traditional alcoholic beverage", "tibetan traditional dairy products",
    "fish meal", "bagged lettuce", "wheat beer", "patagonian merlot wine",
    "apple juice from cider press", "sourdough","shrimp; frozen raw", "palm brown sugar",
    "isolate obtained from campbell's Soup", "brine of stinky tofu", "mead",
    "french dry-type pork sausage", "sichuan pickle vegetables", "emmental",
    "cheese production", "pork dumplings", "home-made koumiss", "korean kefir",
    "chinese egg powder", "stinky tofu","leaf vegetable", "retailed chichen",
    "1 month-old fish sauce mash", "10 weeks old samso 45+ cheese", "10'' wieners",
    "doenjang", "flour", "unsweetened puffed rice cereal", "pinot noir wine",
    "alfalfa sprouts", "ground turkey sausage", "turkey; ground", "sponge cake"
    "jogaejeotgal; a traditional korean", "korean soybean paste",
    "retail chichen", "ground red chili pepper", "makgeollli",
    "raw milk cheese", "maotai daqu", "jeotgal", "zha-chili", "heat treated raw milk",
    "dairy product", "isolate obtained from campbell's soup", "retail pork",
    "traditional homemade dairy product", "gochujang(korean red pepper paste)",
    "tuna steaks; frozen yellowfin", "retail chicken", "retail chichen",
    "chicken carcass", "malted barley", "bee honey", "anaheim pepper",
    "chicken product", "pork steak", "algal food product", "ganjang(korean soy sauce)",
    "bakery environment - concentrated whipped topping",
    "fish part (head) sold as asian bighead carp; muscle",
    "plant derived food stuff; onion; allium cepa", "the brine of stinky tofu",
    "cold-stored modified atmosphere packaged broiler filet strips with the first signs of spoilage"]


livestock_hosts = [
    "sus scrofa", "suxs scrofa domesticus", "bos taurus", "piglet", "pig", "cattle",
    "cows", "broiler chicken", "cow rumen", "bos taurus (bovine)", "horse",
    "sheep placental tissue", "mouse-c57bl/6j", "swine", "chicken", "goose", "duck",
    "bovine","capra hircus", "pet dog", "felis catus", "healthy carrier pig",
    "gallus gallus", "gallus gallus domesticus", "porcine", "cow", "cat", "equus caballus",
    "mare", "calf", "gallus", "canicola", "hen", "meat duck", "diseased pig",
    "equus caballus ferus", "pork", "sus", "bubalus bubalis", "water buffalo",
    "goat", "sheep", "equus ferus caballus", "bos bovis","turkey", "rabbit",
    "bos primigenius taurus", "canine", "laoshan milk goat",
    "canis lupus familiaris","bos taurus coreanae", "yak", "dog", "equine", "ovis aries",
    "young chicken", "gallus gallus domesticus isa15",
    "dairy cow", "poultry", "piglets", "lama glama", "ovis aries (domestic sheep)",
    "duck with tremor",
    "sus domesticus"]


livestock_isolation_sources = [
    "cattle", "chicken", "farm", "turkey feces", "wastewater from pig manure",
    "chicken feces", "nasal cavities of a calf", "bos taurus", "nasal cavities of a pig",
    "pig", "swine nasal swab", "the intestine membrane of a diarrheic piglet",
    "intestine membrane of a diarrheic piglet", "pig manure", "canine oral cavity",
    "isolated from the intestines of sick birds in the farm", "poultry establishment",
    "slaughtered pig", "poultry", "cat", "duck", "broiler", "pig faeces",
    "fecal samples from slaughtered sheep", "broiler chicken","chicken cecal content",
    "dog with mastitis", "goose anus swab", "parent strain cv601 collected from dairy manure",
    "porcine", "cattle hide", "goat; fecal","rumen fluid","pig environment","porcine feces",
    "cow raw milk", "cloacal swab", "skin of a pig", "cattle faeces", "porcine rectal swab",
    "lesion site (lung) of a dead turkey with colibacillosis", "laying hen withcolibacillosis",
    "goose faeces", "cattle slurry", "carcass", "pig fecal", "lamb; fecal",
    "the respiratory tract of a pig with swine respiratory disease","cattle feces",
    "wool from pakistan", "chicken skin", "raw milk", "swine excrement",
    "a nose swab sample of swine origin", "healthy broiler chicken","chicken manure",
    "alpaca; fecal", "horse with strangles", "fresh raw milk", "sus scrofa: cecal content",
    "swine's gut", "swine final chilled carcass", "cow rumen", "cows", "dog",
    "bovine rumen", "cow dung","poultry carcass", "air of cow shed", "pig feces",
    "pig feed from feed plant", "fjerkrae", "liver of poult", "horse feces",
    "tissue and/or biological fluid (swine)", "swine cecum", "raw cow milk",
    "heifer vaginal mucus", "cloaca", "buffalo calf", "aborted piglet fetus",
    "calf", "animal feed", "rectal swab of cattle in slaughterhouse",
    "retail chicken","canine oral plaque","bovine pre-evisceration carcass at",
    "novine (bobby calf; fecal)", "a feces sample of chicken origin",
    "lung of aborted horse fetus", "beef liver", "swine", "chicken intestine",
    "broiler chick cecum","isolated from the intestines of sick", "chicken cecum",
    "bulk pig ears", "cattle slurry", "cow feces", "animal-swine-roaster swine",
    "ncsu equine educational unit", "cowshed of a farm", "a feces sample of swine origin",
    "sheep fecal sample", "digested slurry of dairy manure", "diseased pig",
    "a anal swab sample of swine origin", "ear swab from dog", "goat", "healthy weaning piglets",
    "lesion site (lung) of a dead turkey", "chicken dung", "chicken trachea",
    "pooled sheep faecal samples collected from floor of farm", "poultry litter",
    "pooled pig faecal samples collected from floor of farm", "bull calf intestinal microflora",
    "pooled cattle faecal samples collected from floor of farm",
    "rectal fecal grab samples from a commercial feedlot", "chicken tissues"]


anthropogenic_hosts = ["sewage", "soil (polluted)", "dsmz strain"]


anthropogenic_isolation_sources = [
    "cleaning system aquaculture", "clam larvae aquaculture", "mineral salt medium",
    "combined sewer", "oil", "roundup", "culture mutant", "turbot fish farm",
    "commercial lactic starter", "beet pulp; sugar refinery",
    "high concentration of fluoride", "fish larvae aquaculture", "sewage water",
    "combined sewer effluent", "feed additive imported from china", "wolfram mine tailing",
    "swab from a follow-up assessment of the tap handles and sink edges after a first disinfection attempt",
    "storage tank", "zc4 compost from a compost operation of sao paulo zoo",
    "swab from a hand-washing sink as part of the hospital routine surveillance program",
    "potash salt dump","bakery environment - hallway", "r2a medium", "polystyrene",
    "in the nicotine environment", "outdoor built environment", "farm kitchen",
    "clone isolated from the evolution experiment described in ketola et al. 2013 (doi: 10.1111/evo.12148)",
    "the sample was isolated from cell culture after 9 month of microevolution experiment",
    "salt mine","bioreactor", "sewage & soil", "tung meal", "cooling tower",
    "surface of a polyethylene microplastic particle present in tank 6 of a marine aquarium containing stony-coral fragments and water maintained at 26 degree c",
    "aquarium water", "biofilm sample", "medium", "dairy slurry uk",
    "detritus agregates formed in tank 6 of a marine aquarium containing stony-coral fragments and water maintained at 26 degree c", "duedonoscope instrument", "air in a school dining room",
    "non-filtered water from the water column of tank 6 of a marine aquarium containing stony-coral fragments. water maintained at 26 degree c", "cheddar cheese factory", "mine",
    "anabaena culture", "dairy environment", "genome shuffling", "jugular catheter",
    "telephone of nurse station","shower 3","geothermal power plant","sink aerator",
    "shelf","brewery environment", "lab strain", "coffee cup", "culture collection",
    "air in laboratory", "sp4 medium; universite de rennes", "industrial","knife", "wet market",
    "hospital", "open pond on an algae farm", "carbonate chimney", "concrete",
    "sink handle in icu room in military hospital", "hospital environment",
    "sf9 cell culture media", "laboratory","mold-colonized wall of an indoor",
    "power plant biotrickling filter", "biopesticide dipel df", "biopesticide delfin",
    "biopesticide novodor 3fc", "biopesticide solbac", "biopesticide xentari",
    "alcohol foam dispenser in hospital intensive care unit", "oilfield",
    "bedside light switch in hospital intensive care unit", "1970 vaccine stock",
    "shower 2", "room 7", "research laboratory","mycoplasma culture contaminant",
    "anaerobic digestion reactor", "anthropogenic", "cattle slaughter plant",
    "heated-cooler unit water tank", "beer keg", "woodchip bioreactor",
    "wastewater treatment plant effluent", "ac condensate", "industry",
    "bed sheets","automobile air-conditioning evaporator", "agricultural settling lagoon",
    "biocathode mcl", "envo:01000905", "indoor air", "cultured in erlenmeyer flask",
    "underground farm slurry reservoir", "thermal field", "crude oil",
    "sanger centre via imperial college", "market", "maiket", "cutting board",
    "pulp of a gold-containing sulfide concentrate", "zebrafish tank detritus",
    "polycyclic aromatic hydrocarbon", "agar plate", "brewery", "latex",
    "bakery environment - assembly production room", "hot-gas well; coatings inside tube",
    "anodic biofilm of glucose-fed microbial","milk powder production facility",
    "coke and gas plant treatment facilities", "liquid-air interface biofilm",
    "rubber production plant territory", "laboratory stock","air conditioning system",
    "ensuite 7/8", "river sediment polluted by acid mine drainage",
    "microbial mats from zloty stok gold mine","acid mine decant and tailings from uranium mine",
    "sedimentation pond in a zinc factory","acid mine drainage",
    "orfrc groundwater during biostimulation for uranium bioreduction",
    "sink drain", "chemical manufacturing sites","borehole hdn1; spa", "biofilm boat",
    "e-waste recycling site", "paper mill kaolin", "municipal wwtp",
    "np30 nitrate/perchlorate-reducing mixed", "dust collector of pigpen",
    "acid mining effluent decantation pond", "fish farm",
    "anaerobic digester","cooling tower water", "dairy", "uninoculated hep-2 tissue culture",
    "tattoo ink", "environment of small animal veterinary clinic",
    "arjo bathroom", "paper surface", "composted garbage",
    "floor surface of biological laboratory", "experimental exposure to phage",
    "hot water tap; geest office building", "industrial building air scrubber liquid",
    "iron water pipe", "sterile tools", "dental water pipeline","sungai pinang; penang; malaysia",
    "tar at the shipwreck site of tanker haven", "human septic tank","monitor panel",
    "western north pacific station s1", "ventilator", "coal-cleaning residues",
    "fresh-cut produce processing plant", "burns unit surveillance",
    "chilean kraft-pulp mill effluents", "automobile air conditional evaporator",
    "hospital sink", "oil-immersed sample from guaymas basin", "oil contaninated soil",
    "underground horizons of a flooded mine in russia", "canadian salted buffalo hide",
    "red heat in salted hides (spoiled fish)", "colorado serum co.; anthrax spore vaccine",
    "electroporation of y pestis kim6+ in pcd1ap+", "solar panel array [envo:01000867]",
    "denitrifying; sulfide-oxidizing effluent-treatment plant",
    "high motility on swim plates; streptomycin resistant", "insecticide factory",
    "enrichment cultures from ucc lynda-stan", "from type material",
    "base of single tree in sidewalk along rollins rd. evidence of many dogs; boston; mass",
    "monkey kidney tissue-culture fluids of the fh strain (eaton agent virus)",
    "conjugation assay", "lab", "culture maintained in leon; missing plasmid pscl3",
    "probitic products", "washroom sink in hospital intensive care unit",
    "hygromycin b antibiotic bottle", "paper pulp mill", "residential yard",
    "biofilm reactor", "nutrient broth"]

sediment_isolation_sources = [
    "mud","black mud", "hydrothermal vent area derived sediment", "envo:01001050",
    "marine mud", "muddy water", "pit mud", "tidal marsh", "freshwater mud",
    "mangrove swamp", "sea mud","mud flat", "the pit mud of a chinese liquor", "tidal flat sample",
    "muddy water", "organically rich wetland",
    "mangrove wetland ecosystem", "subseafloor basaltic crust","lake mainaki salt",
    "mud from a salt lagoon"]

fungal_hosts = ["hymenochaete rubiginosa", "himantormia", "veronaeopsis simplex y34",
                "fusarium oxysporum f. sp. cucumerinum", "flammulina filiformis",
                "white-rot fungus phanerochaete chrysosporium", "mushroom"]

fungal_isolation_sources = ["mushroom substrate", "brewery yeast",
                            "white-rot fungus phanerochaete",
                            "hyphae","white-rot fungus phanerochaete chrysosporium"]

human_isolation_sources = [
    "from patient with pneumonia", "blood", "stool","saliva","skin swab", "homo sapiens",
    "wound", "dental abscess", "compound fracture", "dog bite", "oral mucosa",
    "rectal swab", "skin", "fecal", "oral", "nosocomial environment",
    "human tissue biopsy","korean adult feces","groin","infant's throat",
    "sputum (cystic fibrosis patient)", "diarrheal patient", "chest pus",
    "clinical patient","intestine","respiratory tract","human clinical sample",
    "vomit of food poisoning patient", "faecal", "respiratory", "materia alba",
    "human skin swab","vomit from a food poisoning case", "superficial wound",
    "china: the first affiliated hospital of zhejiang university","sputum",
    "human fecal sample","human blood","blood sample","urinary tract infection",
    "bronchial lavage","human feces", "anus swab", "feces; fp", "nasopharyngeal mucus",
    "wound infection","human","stool of a patient presenting with","breast abscess",
    "gastric","infant feces","eye", "semen", "inguinal/rectal swab",
    "human feces (woman; 24 years old)","human stool","urine","faeces","inflamed gingiva",
    "skin abscess", "hepatic abscess", "burn wound", "medical specimen",
    "tibia/osteomyelitis", "patient", "lung sample", "infant faecal",
    "isolated from blood of a patient with","16-year-old boy with bubonic plague",
    "fecal sample","tissue", "ascites", "vaginal secretions", "diarrhea", "genital tract",
    "kidney", "liver", "nose swab", "throat swab", "mouth", "dental plaque", 
    "blood of a hospitalized patient","brain","human feces (woman; 60 years old)",
    "foodborne disease surveillance", "supragingival dental plaque",
    "korean infant feces", "healthy adult male gut","oropharynx","human foot","vein blood",
    "perirectal swab", "rectal", "anal swab", "human skin tissue", "a diarrheal child",
    "blood specimen","rectal fecal grab samples from a","throat", "cystic fibrosis patient",
    "hospitalized patients","from patient with wound infection","human gi tract",
    "human case of meningitis", "supragingival plaque; periodontitis",
    "cancer patient stool", "milker's hand", "hip replacement",
    "biopsy of antral stomach region from european patient with peptic ucler and chronic gastritis disease"]

na_hosts = ["environmental", "jamiecosley", "environment"]

na_isolation_sources = [
    "enviornmental","host's whole body","bacterial consortium", "missing", "unknown",
    "cell culture","afb-diseased colony","environment swab", "physical",
    "biological fluid","not available: to be reported later","environment",
    "r.j.roberts", "whole body", "atcc isolate", "feed", "not available: not collected",
    "isolate obtained from atcc:7955 nca","kimoto", "biofilm", "natural biofilms",
    "obscured","peter dedon mit","swab","dsm:3754","laboratory isolate","tsoundzou", "surface",
    "environmental surface", "ghana", "minnesota", "jiangsu", "not available",
    "usa: wa: seattle","not available: to be reported later","not known","bone powder",
    "laboratory strain derived from cip","local geographical source", "dsmz isolate", "neb433", "bal",
    "air sacs","isolated from commercially available","contaminated wipes","protoplast breeding",
    "iam12617", "singapore", "ncimb", "valley", "conjugation",
    "swab with brown-gray powder", "epidemic",
    "anaerobic environments","environmental","neb269",
    "environmental swabs",
    "integrative microbiology research center; south china agricultural university (scau); tianhe district; guangzhou china",
    "coconi","csf","balf","china: zibo; shandong","wild type strain isolated from a natural source", "beijing",
    "poland; warsaw area","agr",
    "Enviornmental"]

#######################################################################################
#Script Begins Here
#######################################################################################


## CRITICAL TODO: ~650 manual annotations are outdated, so many of them (~450) were marked NA but have been updated to have
## other values. Thus the code does not pick them up properly as it marks some cases as NA when they are not.
## for instance, some new annotations are annotated as "NA", probably due to an upstream bug in which the annotation
## fields changed in the original Genbank/GFF3 annotations. DOUBLE-CHECK ALL THESE CASES!

## For now, use manually-annotated-gbk-annotation-table.csv -- I stripped out the Manual_Annotation
## column to make a "pretend-input" table in order to test this code.

print("Annotation_Accession,host,isolation_source,Annotation")  ## print header

with open("../results/gbk-annotation-table.csv", "r") as annotation_fh:
    for i, line in enumerate(annotation_fh):
        if i == 0: continue ## skip the header.
        line = line.strip()  ## remove trailing newline characters
        fields = line.split(',')  ## split line into an array
        annotation_accession, original_host, original_isolation_source = fields
        ## turn host and isolation source to lower case (NA goes to na!)
        host = original_host.lower()
        isolation_source = original_isolation_source.lower()
        annotation = "blank" ## default value, so we know whether it is set later on.
        
        ## NA
        if host == "na" and isolation_source == "na":
            annotation = "NA" ## the output NA should be upper-case.
        if host == "na" and isolation_source in na_isolation_sources:
            annotation = "NA" ## the output NA should be upper-case.
        elif isolation_source in na_isolation_sources and host == "na":
            annotation = "NA"
        elif host in na_hosts and isolation_source == "na":
            annotation = "NA"
        elif host in na_hosts and isolation_source in na_isolation_sources:
            annotation = "NA"

        ## Terrestrial
        elif host == "na" and isolation_source in terrestrial_isolation_sources:
            annotation = "Terrestrial"
        elif "rifle well" in isolation_source:
            annotation = "Terrestrial"
        elif "soda lake" in isolation_source:
            annotation = "Terrestrial"
        elif "hypersaline" in isolation_source:
            annotation = "Terrestrial"
        elif "hot spring" in isolation_source and "soil" not in isolation_source:
            annotation = "Terrestrial"
        elif host == "air":
            annotation = "Terrestrial"

        ## Soil
        elif "soil" in isolation_source and "contamina" not in isolation_source and "pollute" not in isolation_source and "contaninated" not in isolation_source and "sewage" not in isolation_source:
            annotation = "Soil"
        elif host == "na" and isolation_source in soil_isolation_sources:
            annotation = "Soil"
        elif "soil" in host and "contamina" not in host and "pollute" not in host:
            annotation = "Soil"
        elif "rhizosphere" in isolation_source or "rhizosphere" in host:
            annotation = "Soil"

        ## Agriculture
        elif host in agri_hosts and host not in plant_hosts and isolation_source not in food_isolation_sources:
            annotation = "Agriculture"
        elif isolation_source in agri_isolation_sources:
            annotation = "Agriculture"
        elif "sativa" in host:
            annotation = "Agriculture"
        elif "citrus" in host:
            annotation = "Agriculture"
        elif "solanum" in host:
            annotation = "Agriculture"
        elif "tomato" in host:
            annotation = "Agriculture"

        ## Plant-host
        elif host in plant_hosts:
            annotation = "Plant-host"
        elif isolation_source in plant_isolation_sources:
            annotation = "Plant-host"
            
        # Marine
        elif host in marine_hosts:
            annotation = "Marine"
        elif isolation_source in marine_isolation_sources:
            annotation = "Marine"
        elif host == "na" and " sea " in isolation_source and "sediment" not in isolation_source and "mud" not in isolation_source and "yogurt" not in isolation_source:
            annotation = "Marine"
        elif "seawater" in isolation_source and host == "na" and "aquaculture" not in isolation_source:
            annotation = "Marine"
        elif "hydrothermal" in isolation_source and host == "na" and "sediment" not in isolation_source:
            annotation = "Marine"
        elif ("marine" in isolation_source and "sediment" not in isolation_source and "mud" not in isolation_source
              and host not in animal_hosts and isolation_source not in terrestrial_isolation_sources
              and "sponge" not in isolation_source and "aquarium" not in isolation_source):
            annotation = "Marine"

        ## Freshwater
        elif host in freshwater_hosts:
            annotation = "Freshwater"
        elif host == "rhodobacter sphaeroides 2.4.1":
            annotation = "Freshwater"
        elif ("water" in isolation_source and isolation_source not in anthropogenic_isolation_sources
              and isolation_source not in sediment_isolation_sources
              and "sea" not in isolation_source and "waste" not in isolation_source and "oil" not in isolation_source and "aquarium" not in isolation_source
              and host == "na" and "ocean" not in isolation_source and "sediment" not in isolation_source
              and isolation_source not in food_isolation_sources and "sal" not in isolation_source and "kimchi" not in isolation_source):
            annotation = "Freshwater"
        elif "spring" in isolation_source and "oil" not in isolation_source and "sediment" not in isolation_source:
            annotation = "Freshwater"
        elif ("lake" in isolation_source and "sediment" not in isolation_source and "sal" not in isolation_source
              and "soda" not in isolation_source and "oil" not in isolation_source and "silt" not in isolation_source):
            annotation = "Freshwater"
        elif isolation_source in freshwater_isolation_sources:
            annotation = "Freshwater"

        # Food
        elif host in food_hosts:
            annotation = "Food"
        elif "pickle" in host:
            annotation = "Food"
        elif isolation_source in food_isolation_sources:
            annotation = "Food"
        elif (host == "na" and "milk" in isolation_source and "hand" not in isolation_source
              and "raw" not in isolation_source and "facility" not in isolation_source):
            annotation = "Food"
        elif host == "na" and "food" in isolation_source and "vomit" not in isolation_source and "sludge" not in isolation_source:
            annotation = "Food"
        elif host == "na" and "beef" in isolation_source:
            annotation = "Food"
        elif "cheese" in isolation_source:
            annotation = "Food"
        elif "pickle" in isolation_source:
            annotation = "Food"
        elif "kimchi" in isolation_source:
            annotation = "Food"
        elif "meat" in isolation_source and "facility" not in isolation_source:
            annotation = "Food"
        elif "ferment" in isolation_source:
            annotation = "Food"
        elif host == "na" and "chicken" in isolation_source and isolation_source not in livestock_isolation_sources:
            annotation = "Food"
        elif host == "na" and "yogurt" in isolation_source:
            annotation = "Food"

        ## Livestock
        elif host in livestock_hosts and isolation_source not in food_isolation_sources:
            annotation = "Livestock"
        elif isolation_source in livestock_isolation_sources and host not in animal_hosts:
            annotation = "Livestock"
        elif "bovine" in isolation_source:
            annotation = "Livestock"
        elif "scrofa" in host:
            annotation = "Livestock"
        elif "placental tissue" in isolation_source:
            annotation = "Livestock"
        elif isolation_source == "turkey" and host == "na":
            annotation = "Livestock"

        ## Animal-host
        elif host in animal_hosts:
            annotation = "Animal-host"
        elif isolation_source in animal_isolation_sources:
            annotation = "Animal-host"

        ## Anthropogenic-environment
        elif host in anthropogenic_hosts:
            annotation = "Anthropogenic-environment"
        elif host == "na" and isolation_source in anthropogenic_isolation_sources:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "aquaculture" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "waste" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "sewage" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "sludge" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "contamina" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "digester" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif "oil" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif "cell line" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "industr" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "facility" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "engineer" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "reactor" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "laboratory" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "experiment" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "biopesticide" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif "polluted soil" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "gram-positive bacteria" or host == "bacteria" or host == "sludge":
            annotation = "Anthropogenic-environment"
        elif "probiotic" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host ==  "sofa":
            annotation = "Anthropogenic-environment"
            
        ## Sediment
        elif "sediment" in isolation_source and "contamina" not in isolation_source and "pollute" not in isolation_source:
            annotation = "Sediment"
        elif host == "na" and isolation_source in sediment_isolation_sources:
            annotation = "Sediment"
        elif host == "marine sediment":
            annotation = "Sediment"
        elif "mud" in isolation_source:
            annotation = "Sediment"
        elif "silt" in isolation_source:
            annotation = "Sediment"
        elif "wetland" in isolation_source:
            annotation = "Sediment"
        elif "estuary" in isolation_source:
            annotation = "Sediment"
            
        ## Fungal-host
        elif host in fungal_hosts or isolation_source in fungal_isolation_sources:
            annotation = "Fungal-host"

        ## Human-host
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and "homo" in host:
            annotation = "Human-host"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and "human" in host:
            annotation = "Human-host"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and host == "human":
            annotation = "Human-host"
        elif host not in livestock_hosts and host not in animal_hosts and isolation_source in human_isolation_sources:
            annotation = "Human-host"
        elif host == "na" and isolation_source == "feces":
            annotation = "Human-host"
        elif host == "na" and "blood" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and "patient" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and "human" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and "stool" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and "clinical" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and isolation_source == "lung":
            annotation = "Human-host"

            
        print(','.join([annotation_accession, original_host, original_isolation_source, annotation]))

