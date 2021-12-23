#!/usr/bin/env python

'''
annotate-ecological-category.py by Vincent Huang and Rohan Maddamsetti.

This script annotates the ecology for each analyzed genome, based on the metadata in the
'host' and 'isolation_source' tags in its genome file from Genbank.

Some definitions for annotation: soil is everything from dirt, whether an agricultural field or not.
agriculture, on the other hand, refers strictly to crop plants-- soil is not included.

In addition, "rhizosphere" always maps to the Soil annotation.

Usage: python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv
'''

########################################################################################################################
## Keywords for Annotating Ecological Categories.
########################################################################################################################

animal_hosts = [
    "donacia thalassina", "ixodes persulcatus", "tegillarca granosa", "wild bird",
    "edessa sp.","nematode", "seriola lalandi", "marmota himalayana", "egret",
    "tragelaphus strepsiceros", "apis mellifera", "botryllus sp.", "oncorhynchus kisutch",
    "meleagris gallopavo", "blue sheep", "shellfish", "insect", "chelymorpha alternans",
    "neophocaena asiaeorientalis", "penaeus japonicus", "alligator", "macroplea mutica",
    "tibetan antelope", "penaeus vannamei (whiteleg shrimp)", "larus sp.", "protozoa",
    "donacia provostii","donacia marginata", "ictalurus punctatus", "amblyomma cajennense",
    "argas persicus","mink", "cistudinella sp.", "ischnocodia annulus", "bigeye tuna",
    "drosophila melanogaster oregon-r-modencode", "sablefish", "oncorhynchus mykiss", "atlantic salmon",
    "tilapia fish", "migratory bird", "south china tiger", "crocodile lizard",
    "arvelius albopunctatus", "crassostrea gigas","phoca largha","amblyomma variegatum (cattle tick)",
    "equus kiang", "phractocephalus hemioliopterus", "black-collared starling", "marmot",
    "cryptocercus clevelandi", "wild chukar", "alvinocaris longirostris",
    "penaeus vannamei","donacia versicolorea","nipponaphis monzeni","salmo salar","amazona sp.",
    "plateau pika", "silkworm", "gentoo penguin", "discomorpha sp.",
    "chroicocephalus novaehollandiae (australian silver gull chick)", "trachinotus ovatus", "shrimp",
    "crow", "cassida sp.", "pelodiscus sinensis", "hirudo verbana",
    "galaxea fascicularis (stony coral)","bivalve mollusk","ondara zibethicus","aphis glycines",
    "mus musculus", "sebastes schlegeli", "mus musculus subsp. domesticus",
    "donacia semicuprea","tuberolachnus salignus", "hippopotamus amphibius", "paralichthys olivaceus",
    "bos mutus", "pheasant", "charidotella sexpunctata",
    "neoaliturus tenellus", "crucian carp", "ictalurus furcatus", "mercenaria mercenaria",
    "physopelta gutta", "agroiconota sp.", "pentalonia nigronervosa",
    "apterostigma dentigerum","urocitellus undulatus", "mussels", "artemisaphis artemisicola",
    "aphis nerii","eel", "ixodes ricinus (tick)", "acromis sparsa",
    "macroplea appendiculata","donacia clavipes","connochaetes taurinus","perinereis linea",
    "trichonympha agilis", "wild boar", "pseudotrichonympha sp.",
    "drosophila melanogaster","dicentrarchus labrax (sea bass)", "rainbow trout", "donacia dentata",
    "honey bee", "channa argus (snakehead fish)",
    "ixodes scapularis","gut of wasp", "american cockroach", "trypoxylus dichotomus", "oyster larvae",
    "parrot", "paper wasp", "marmota", "scophthalmus maximus",
    "oreochromis niloticus", "ornithodoros hermsi", "oecophylla smaragdina", "euschistus servus",
    "ixodes ricinus", "cassida viridis", "cassida vibex",
    "phakellia ventilabrum", "phoca vitulina", "megacopta punctatissima", "invertebrates",
    "acyrthosiphon pisum", "donacia vulgaris","mouse", "mya arenaria oonogai makiyama",
    "nasonia vitripennis","cryptocercus punctulatus",
    "schizaphis graminum biotype i", "oreochromis", "cage-cultured red drum", "coral",
    "ornithodoros sonrai", "fish", "amblyomma longirostre",
    "gorilla", "hydra vulgaris aep", "amblyomma neumanni","plateumaris rustica","oryctes gigas",
    "crab","honeybees", "aratinga solstitialis",
    "syngnathus typhle","nezara viridula (cotton pathogen vector; southern green stink bug)",
    "litopenaeus vannamei","naegleria",
    "deinagkistrodon acutus", "aphis helianthi", "anguilla anguilla", "halyomorpha halys",
    "antho dichotoma", "grass carp",
    "pthirus gorillae","marmota baibacina","ruddy shelduck", "marmota sibirica",
    "protaetia brevitarsis seulensis",
    "bactrocera oleae","rattus", "hyadaphis tataricae", "myzus persicae", "alces alces",
    "flatfish","toothfish", "bankia setacea",
    "acyrthosiphon kondoi", "sebastes schlegelii", "penaeus vannamei (shrimp)", "ovis",
    "macrosiphum gaurae",
    "turbot", "pseudotrichonympha grassii (protist) in the gut of the termite coptotermes formosanus",
    "cerambycidae sp.", "pediculus humanus corporis",
    "kangaroo", "cryptopsaras couesii","anas platyrhynchos", "penaeus setiferus",
    "apis mellifera mellifera",
    "camponotus chromaiodes", "caenorhabditis elegans","murgantia histrionica",
    "ixodes pacificus","donacia proxima",
    "parakeet", "geronticus eremita","melanaphis sacchari","dendrolimus ibiricus",
    "drosophila melanogaster oregon", "allacta bimaculata",
    "myzus persicae (green peach aphid)", "macrotermes barneyi", "glossina morsitans morsitans",
    "nauphoeta cinerea", "delphinapterus leucas (beluga whale)",
    "oryctolagus cuniculus","corvus brachyrhynchos","drosophila melanogaster oregon-r modencode",
    "apterostigma",
    "odocoileus virginianus","peromyscus leucopus", "cinara cedri", "brachycaudus cardui",
    "seabass", "rabbit",
    "plateumaris consimilis", "acanthamoeba polyphaga hn-3", "perameles bougainville",
    "plateumaris braccata","rat",
    "hyalomma aegyptium","zebra", "apsterostigma", "baizongia pistaciae", "pediculus schaeffi",
    "sitobion avenae",
    "otospermophilus beecheyi (california ground squirrel)","apodemus agrarius", "donacia fulgens",
    "bat", "skate",
    "paralichthys olivaceus (flounder)","dermacentor andersoni","donacia tomentosa","lacertilia",
    "papio papio",
    "anomalops katoptron", "gull", "japanese eel", "salvelinus fontinalis","blattella germanica",
    "pteropus livingstonii", "pteropus poliocephalus",
    "moschus berezovskii", "macaca silenus", "panulirus ornatus","pavona duerdeni",
    "odontobutis platycephala",
    "ostrea edulis (flat oyster)", "beluga whale", "shinkaia crosnieri",
    "lutjanus guttatus (rose snapper)",
    "vultur gryphus","penaeus monodon","tanakia koreensis", "dicentrarchus labrax",
    "phyllophaga sp.","bombyx mori",
    "zophobas atratus", "pigeon", "turtle", "penaeus (litopenaeus) vannamei (whiteleg shrimp)",
    "donacia cinerea", "fulmars", "brevicoryne brassicae", "acanthamoeba", "blaberus giganteus",
    "arctic char", "microlophium carnosum",
    "macaca mulatta","thelaxes californica","bactrocera dorsalis","apostichopus japonicus",
    "ovis aries",
    "pthirus pubis","manis javanica (pangolin)", "sipalinus gigas", "geodia barretti",
    "macrosiphoniella sanborni",
    "crassostrea gigas (pacific oyster)", "crassostrea virginica", "ellychnia corrusca",
    "danaus plexippus",
    "glossina brevipalpis","marine sponge lissodendoryx isodictyalis in the bahamas",
    "muscaphis stroyani","aphis urticata",
    "catfish", "morone chrysops x morone saxatilis", "pediculus humanus", "rattus norvegicus",
    "plateumaris sericea",
    "bothriocroton concolor", "hyperomyzus lactucae","trimyema compressum", "macaca fascicularis",
    "diplonemea",
    "neohaemonia nigricornis", "haemaphysalis juxtakochi", "columba livia", "diaphorina citri",
    "prawn","deer",
    "pseudorca crassidens","acyrthosiphon lactucae","aphis nasturtii","giant panda",
    "cryptocercus kyebangensis",
    "japanese rhinoceros beetle larva","hydrophilus acuminatus", "discus", "hawk",
    "donacia sparganii", "oyster", "caprine",
    "acinonyx jubatus","chaeturichthys stigmatias", "donacia cincticornis", "anas strepera",
    "pan troglodytes", "paguma larvata",
    "aphis craccivora (cowpea aphid)", "calanoid copepod", "atlantic white-sided dolphin",
    "ctenocephalides felis", "parachirida sp.",
    "pyropia yezoensis conchocelis", "pan paniscus", "lipaphis pseudobrassicae",
    "donacia piscatrix", "diaphorina cf. continua",
    "aphis craccivora","donacia crassipes", "donacia bicoloricornis", "donacia simplex",
    "draeculacephala minerva",
    "anodonta arcaeformis", "mastotermes darwiniensis", "toxic alexandrium minutum",
    "coreoleuciscus splendidus", "avian",
    "platycercus elegans", "solea senegalensis","sanzinia madagascariensis volontany",
    "nezara viridula", "lion-tailed macaques",
    "plateumaris pusilla","gilthead seabream", "ixodes pacificus (western blackleg tick)",
    "mus musculus c57bl/6j", "tick",
    "coregonus clupeaformis (lake whitefish)", "wax moth", "galleria mellonella",
    "cyclopterus lumpus","cockatiel",
    "blatta orientalis", "seriola dumerili", "white stork", "canis latrans","drosophila neotestacea",
    "shrimps", "bird of prey",
    "c57bl/6ntac mice", "phascolarctos cinereus", "neopsylla setosa", "haliclona simulans",
    "ardea cinerea",
    "rhopalosiphum padi", "silurus asotus", "uroleucon ambrosiae","mustela putorius furo",
    "melanocetus johnsonii",
    "ailuropoda melanoleuca","spermophilus sp.", "squirrel", "dendrolimus sibiricus", "cormorant"]

animal_isolation_sources = [
    "animal", "wild rat", "coho salmon", "fruit fly","diaphorina citri", "south china tiger",
    "iguana", "bat", "wild yak feces", "white-crowned sparrow", "african elephants",
    "efb-infected honey bee colony", "turban shell", "mussels", "rainbow trout",
    "foulbrood of honeybees", "turbot", "red fox", "panda",
    "marine sponge", "chamaeleonidae","ixodes spinipalpis", "barnacle at wood pile-on",
    "atlantic salmon", "manis javanica", "a feces sample of migratory birds origin",
    "tilapia", "wild pig; fecal", "shellfish","diarrheal snake diarrheal snake in hunan",
    "diseased shrimp", "wasp honeycombs", "eclectus roratus feces",
    "insect larvae", "tissue; animal","yellowtail", "snake", "hydrophilus acuminatus",
    "crassostrea gigas","pigeon", "gill tissue of bathymodiolus japonicus",
    "meleagris gallopavo","sea bass","mouse gut","avian","sockeye salmon",
    "diseased labeo rohita fish","mealworm", "feces; gull (larus spp.)",
    "musca domestica","atlantic cod","siniperca scherzeri", "red-breasted parakeet",
    "midgut crypts of stink bug togo hemipterus","konosirus punctatus",
    "crushed cell of callyspongia sp.", "tuna", "housefly", "bison", "coral","catfish",
    "mactra veneriformis", "black-collared starling", "white-lipped deer",
    "shrimp", "black-headed gull", "psittacus erithacus feces",
    "separated from the corpses of silkworms that had died due to bb natural infection in daiyue district; taian city; shandong province; china.", "animal hide","tuna scrape; yellowfin",
    "tuna scrape; frozen yellowfin nakaoachi",
    "aphis craccivora (cowpea aphid) on robinia pseudoacacia (locust)", "swiss alpine ibex feces",
    "bearded dragon", "wild bird", "rattus rattus",
    "neophocaena phocaenoides", "loach","unidentified actinians", "fish",
    "waxworms gut", "flea", "intestinal contents of termite nasutitermes nigriceps",
    "muskrat during outbreak"]

plant_hosts = [
    "suaeda salsa","cladonia borealis","medicago truncatula", "ulva pertusa (algae)",
    "allamanda cathartica","astragalus pelecinus","lebeckia ambigua","lotus corniculatus",
    "trifolium repens", "curcuma aromatica", "musa balbisiana cultivar kepok",
    "datisca glomerata","beta vulgaris","robinia pseudoacacia","ficus religiosa l.",
    "sporobolus anglicus", "casuarina equisetifolia",
    "leontopodium alpinum","tanacetum vulgare","catharanthus roseus","stereocaulon sp.",
    "oxytropis pumilio", "malus sylvestris", "angelica sinensis dlies",
    "agave americana l.","stachytarpheta glabra","elymus tsukushiensis","eucalyptus",
    "androsace koso-poljanskii", "bromus inermis",
    "trifolium pratense","lemna minor","mimosa affinis","salsola stocksii", "mimosa",
    "biserrula pelecinus l.","kaempferia sp. rhizome","masson pine",
    "nerium oleander","biserrula pelecinus","lichen","miscanthus giganteus","mimosa scabrella",
    "rosa sp.", "malus sieversii", "onobrychis viciifolia",
    "flower","mimosa flocculosa","trifolium uniflorum","eucalypti of eucalyptus","wild cotton",
    "urochloa reptans", "cotinus coggygria", "kandelia candel (mangrove)",
    "trifolium spumosum l. (annual mediterranean clovers)","microcystis culture",
    "medicago orbicularis","arabidopsis", "vachellia farnesiana",
    "physcomitrella patens","cryptomeria fortunei","festuca arundinacea","oxytropis triphylla",
    "usnea","rubus sp.", "melilotus officinalis",
    "achillea ptarmica l.","muscari","lespedeza cuneata","nerium oleander",
    "tanacetum vulgare","eucommia ulmoides",
    "leiosporoceros dussii","oxytropis triphylla","vavilovia formosa", "oxytropis kamtschatica",
    "zantedeschia aethiopica","rosa", "prunus cerasifera",
    "leontopodium nivale","nusuttodinium aeruginosum","prosopis cineraria",
    "oryza glumipatula","acacia farnesiana",
    "salvia splendens","panax ginseng", "vicia alpestris","phaeoceros",
    "commiphora wightii", "alhagi sparsifolia shap.", "acacia farnesiana",
    "lobaria pulmonaria thallus", "lotus",
    "hyacinthus orientalis","sporobolus anglicus","malus prunifolia (crab apple)","tephrosia apollinea",
    "aeschynomene fluminensis","vavilovia formosa", "prunus sp.",
    "oenothera speciosa", "acaciella angustissima", "brassica rapa subsp. chinensis",
    "catharanthus roseus","dracaena sanderiana","parthenium argentatum gray (guayule shrubs)",
    "medicago arborea", "medicago lupulina l.",
    "euonymus sp.","eucalyptus grandis","broussonetia papyrifera","lichen stereocaulon sp.",
    "lotus sp.","weed", "thlaspi arvense", "grass",
    "blasia pusilla","calliandra grandiflora","brittle root","digitaria eriantha","vavilovia formosa",
    "plant", "sesame", "salix sp.", "hibiscus"]

plant_isolation_sources = [
    "flower of forsythia koreana", "algal phycosphere", "lotus corniculatus nodule", "malus sylvestris",
    "flower of rhododendron schlippenbachii", "root of codonopsis pilosula", "root nodules",
    "chrysochromulina tobin phycosphere", "algae","artificially infected catharanthus", "leaf spots",
    "peltigera membranacea thallus","angelica gigas nakai root surface",
    "galega orientalis root nodule", "seaweed", "xylem",
    "plant","broussonetia papyrifera","lichen", "lowbush blueberries",
    "seaweed; enteromorpha linza (coastal marine sediments); aburatsubo inlet; australia",
    "flower","leaf of fig tree", "grass from the park in a vetenary clinic",
    "plant roots","galega officinalis root nodule","lichen thallus","ginseng","gladiolus",
    "sesbania spp. root nodule","plant leaf","origanum marjorana","malus sp.","algae",
    "flower of shittah tree",
    "inner tissues of halophyte limonium sinense (girard) kuntze","leaf",
    "root nodule of sesbania cannabina",
    "flower of rhododendron sclippenbachii", "plant xylem"]

terrestrial_isolation_sources = [
    "iron hydroxide deposits", "solar salt", "air", "elton hypersaline lake",
    "brazilian saline-alkaline lake", "desert", "saltpan",
    "coal bed","moderate hot spring","salt marsh","alkaline pool submerged anode electrode",
    "kulunda steppe hypersaline lake", "hypersaline environment",
    "chronically low temperature and dry polar region",
    "solar saltern", "dust", "beach sand", "salar de atacama; atacama desert","weathered rock sample",
    "moist arsenopyrite (feass)-containing rock taken from a mine tunnel approximately 300 m below the ground in the granites gold mine", "microbial mat",
    "inside the caves of drach", "non-purified solar salt","gold-copper mine",
    "xinjiang aibi salt lake",
    "shar-burdin hypersaline soda lake", "marine solar saltern brine","soda lake magadi",
    "rock (ENVO:00001995)",
    "lava", "brine", "saltern", "baengnokdam summit crater area; mt. halla",
    "taibei marine solar saltern near lianyungang city","orthoquartzite cave surface","cave",
    "the bange salt-alkaline lake in tibet",
    "geothermal isolate", "salt crust", "coal spoil heap",
    "flat; laminated microbial mat in a salt marsh",
"rock (envo:00001995)", "saline saltern", "hypersaline lake",
    "thiodendron' bacterial sulfur mat from mineral sulfide spring"]

soil_isolation_sources = [
    "agricultural field", "tundra wetland", "permafrost", "compost", "mixed sand sample",
    "sphagnum peat from the bog obukhovskoe (acidic wetland)","paddy field; sungai manik; malaysia",
    "rice fields", "sphagnum peat", "soil around hot spring", "dirt",
    "sphagnum bog", "fertilizer", "ancient permafrost from mammoth","montane grasslands",
    "solar salt farm", "enriched culture of compost",
    "siberian permafrost","particulate matter", "pine forest", "mushroom compost", "manure compost",
    "composted cattle manure", "agricultural waste material","long-term organic manure fertilized",
    "Microscale soil grain", "ermafrost region of qilian mountains"]

agri_hosts = [
    "coffea", "grape", "citrus spp.", "prunus avium", "capsicum sp.","phaseolus vulgaris",
    "sugarcane", "lolium perenne", "zingiber officinale",
    "phaseolus", "punica granatum", "pyrus pyrifolia var. culta","alfalfa","cucumis sativus",
    "cucurbita maxima", "lentil", "lathyrus sativus",
    "solanum melongena", "vitis vinifera","vitis vinifera l. cv. seto giants", "grapevine",
    "sakura tree", "rice", "japanese radish",
    "chinese cabbage","kiwifruit", "pisum sativum l. (pea)", "common bean", "melon", "pepper",
    "glycine max", "glycine", "manihot esculenta",
    "actinidia deliciosa", "sugarcane", "actinidia chinensis", "lettuce", "olea europaea",
    "grapefruit", "brassica juncea var. foliosa", "apple",
    "brassica oleracea var. botrytis", "camellia sinensis",
    "brassica rapa subsp. pekinensis",
    "triticum aestivum", "areca catechu", "date palm", "zea mays l.", "ginger","plum", "cowpea",
    "potato", "wheat",
    "saccharum officinarum", "cotton", "soybean","mandarin orange", "fava bean", "capsicum annuum",
    "vitis vinifera l. cv. aurora black","ziziphus mauritiana lam", "pyrus pyrifolia", "strawberry",
    "sesame seedling", "tobacco","peanut", "carrot", "apple tree", "banana",
    "soybean (glycine max (l.) merrill)",
    "panicum miliaceum", "solanum lycopersicum", "plantain", "pear", "winter wheat",
    "pyrus communis 'williams'",
    "actinidia", "zea mays", "turfgrass", "prunus dulcis","glycine soja", "potato",
    "pepper plant", "sweet orange",
    "cicer arietinum","pisum sativum", "actinidia deliciosa 'hayward'", "arachis hypogaea",
    "forage rape", "pear",
    "prunus cerasus (sour cherry)", "prunus dulcis", "sorghum bicolor", "phaseolus sp.",
    "allium cepa",
    "vitis vinifera cv. 'izsaki sarfeher'","brassica rapa var. laciniifolia subvar. oblanceolata",
    "allium cepa l.",
    "gossypium hirsutum","solanum lycopersicoides", "vicia faba", "eggplant", "pear tree",
    "juglans regia",
    "glycine max cv. jinju1", "lycopersicon esculentum", "camellia oleifera", "allium cepa (onion)",
    "musa sp.",
    "triticum aestivum l.", "coffee plant", "pineapple", "citrus sinensis", "gossypium sp.",
    "musa spp."]

agri_isolation_sources = [
    "ucb-1 pistachio rootstock", "pineapple", "plants with moko disease","mangoes", "chives",
    "raw almonds", "chilli", "eggplant", "cucumber", "rice leaves",
    "isolated from cucumber", "cabbage","almond drupe", "mexican lime leaf", "apple twig",
    "phaseolus vulgaris root nodule",
    "healthy tomato plant","almond kernel (raw; variety carmel)", "the root of rice", "soybeans",
    "wheat germ", "raw peanuts", "rot potato tubers",
    "sugarcane root", "originally isolated from olive trees","isolated from the tobacco substrate",
    "orange tree", "rice shoot",
    "hydroponic pots with potatoes", "ragi", "canola roots","leaf from rome apple cultivar",
    "cilantro","rice root",
    "orange","maize stem", "strawberry leaf tissue", "rye silage", "fodder","stable grass silage",
    "beans","carrot", "potato",
    "corn", "japanese pear","strawberry leaf tissue", "wheat root", "ogi (red sorghum)",
    "campbell early grape",
    "oryza sativa", "nodules from common bean", "pisum sativum root-nodule",
    "naturally-infected soybean tissue", "soybean",
    "naturally-infected soybean leaf tissue", "apple tree", "banana", "cabbage seeds",
    "citrus orchard", "tomato",
    "white nectarines fire pearl variety"]

marine_hosts = ["argopecten purpuratus", "seawater", "pyropia tenera",
                "trichodesmium erythraeum ims101", "asterionellopsis glacialis strain a3",
        "rhodosorus marinus", "red alga", "pyropia"]

marine_isolation_sources = [
    "salt lake", "hydrothermal vent area derived", "deep sea",
    "water from the baltic sea", "envo:00000227", "sea", "envo:00002149",
    "sea water", "deep-sea hydrothermal vent", "sea ice", "oceanic water",
    "cold seep", "envo:01000301",
    "ocean", "isolated from the deepest ocean",
    "chrysochromulina tobin phycosphere", "gulf of finland", "cyanobacterial bloom",
    "glass slide place on reef flat in natural environment",
    "surface water of the southern north sea","antarctic iceberg",
    "hydrothermal vent", "acidic salty water",
    "shallow tropical waters; normally from coral reef substrate","open ocean water",
    "sulfidic waters (60 m) from the peruvian upwelling region","shellfish hatchery",
    "bay of bengal", "pacific plankton",
    "culture of roseovarius sp. tm1035", "coastal","runway 10 reef (10-12 m)",
    "german wadden sea",
    "microbial mat material from brackish estuary", "hydrothermal precipitates",
    "shallow sea; symbiosis with algae",
    "beach sand","estuarine water [envo:01000301]", "150 km offshore",
    "cyanobacterial aggregates","synechocystis sp. gt-l",
    "seawater", "shallow-sea hydrothermal system"]

freshwater_isolation_sources = [
    "creek","pond","tianshan glacier",
    "glomma river","snow","microcystis culture","glacial stream",
    "envo:00002011", "stream biofilm",
    "euglena gracilis from city ponds", "glacier from lahual spiti valley","river",
    "pool at botanical garden",
    "ayakekumsalt lake", "envo:00002006"]

food_hosts = ["kimchi cabbage", "mung bean", "kombucha scoby", "seujeot", "pickle", "pickled cabbage",
              "beef", "koumiss", "dairy cattle", "nostoc flagelliforme", "sichuan pickle"]

food_isolation_sources = [
    "saeng-gimol meju",
    "futsai", "koumiss", "soybean paste (chonggugjang)", "mashed potatoes", "fermented food",
    "kombucha scoby", "kefir", "cheese product", "instant pork", "pork &cabbage dumplings",
    "isolated from commercially available caspian sea yogurt",
    "makgeolli (korean traditional alcoholic beverage)", "cream products",
    "light wheat beer", "sprouts", "3-day-old traditional semihard cheese","pork",
    "moto starter of sake", "minced pork", "shrimp paste",
    "isolated from salted brown alga laminaria", "shima dofu (okinawan-style tofu)",
    "salted brown alga laminaria", "stinky xiancaigeng", "sake(hatsuzoe)",
    "blue cheese in wax", "tibicos symbiotic community", "pickled cabbage",
    "ready to eat mixed salad leaves (obtained from discount store)",
    "degraded sugar thick juice","fresh alfalfa sprouts", "scallop", "chinese sausages",
    "ground turkey", "gimbap", "product-eggs-raw-whites",
    "frozen peas", "palm wine", "natural whey culture from gruyere cheese","homemade koumiss",
    "salad", "rotten eggs", "steamed conch", "rice wine rice syrup",
    "water kefir", "palm sap","peanut butter", "ground turkey", "sick cider", "onion",
    "austrian hard cheese rind", "litopenaeus vannamei purchased at the supermarket",
    "chinese pickle", "cheese", "raw mutton","nuruk; korean traditional beverage starter",
    "crab marinated in soy sauce",
    "honey", "ham", "bobby veal steak", "10' wieners", "natto", "yogurt","gochujang", "cheonggukjang",
    "traditional greek kasseri cheese","pork & cabbage dumplings", "lettuce", "new zealand cheese",
    "beer", "tibetan kefir", "carrot juice", "pickle", "butter starter", "product-eggs-raw-whole",
    "raw sausage", "beer contaminant", "scoby from kombucha tea", "calf liver","mustard pickles",
    "smoked fish",
    "doubanjiang", "retail turkey","fresh produce (lettuce; lollo bionda)", "blue berry",
    "apple cider","fried eel",
    "mac and cheese", "shelled pistachios", "marinated fish product","10 weeks old 45+ samso cheese",
    "salami",
    "mixed salads", "tibet kefir","sugar thick juice", "traditional dairy products", "ice cream",
    "water containing garlic and cabbage", "doenjang(soybean paste)","green chile ingredient",
    "ground beef",
    "koumiss (fermented mare's milk)", "pilsner beer","blown cheese", "pickled green chili peppers",
    "wheat sourdough", "cherry","pork barbeque", "vinegar", "dried tofu", "sauerkraut", "diced lamb",
    "salt",
    "korean traditional alcoholic beverage", "tibetan traditional dairy products","fish meal",
    "bagged lettuce",
    "apple juice from cider press", "sourdough","shrimp; frozen raw", "palm brown sugar",
    "wheat beer",
    "isolate obtained from campbell's Soup", "brine of stinky tofu", "mead",
    "french dry-type pork sausage", "sichuan pickle vegetables", "emmental","cheese production",
    "pork dumplings",
    "chinese egg powder", "stinky tofu","leaf vegetable", "retailed chichen", "home-made koumiss",
    "korean kefir",
    "1 month-old fish sauce mash", "10 weeks old samso 45+ cheese", "10'' wieners", "doenjang",
    "alfalfa sprouts", "ground turkey sausage",
    "jogaejeotgal; a traditional korean", "korean soybean paste", "retail chichen",
    "raw milk cheese", "maotai daqu", "jeotgal",
    "dairy product", "isolate obtained from campbell's soup", "retail pork",
    "traditional homemade dairy product",
    "tuna steaks; frozen yellowfin", "retail chicken", "retail chichen", "chicken carcass",
    "chicken product", "pork steak",
    "bakery environment - concentrated whipped topping",
    "plant derived food stuff; onion; allium cepa",
    "cold-stored modified atmosphere packaged broiler filet strips with the first signs of spoilage"]


livestock_hosts = [
    "sus scrofa", "suxs scrofa domesticus", "bos taurus", "piglet", "pig", "cattle",
    "cows", "broiler chicken", "cow rumen", "bos taurus (bovine)", "horse",
    "sheep placental tissue", "mouse-c57bl/6j", "swine", "chicken", "goose", "duck",
    "bovine","capra hircus", "pet dog", "felis catus",
    "gallus gallus", "gallus gallus domesticus", "porcine", "cow", "cat", "equus caballus",
    "mare", "calf", "gallus", "canicola",
    "equus caballus ferus", "pork", "sus",
    "goat", "sheep", "equus ferus caballus", "bos bovis","turkey", "rabbit",
    "bos primigenius taurus", "canine",
    "canis lupus familiaris","bos taurus coreanae", "yak", "dog", "equine", "ovis aries",
    "young chicken", "gallus gallus domesticus isa15",
    "dairy cow", "poultry", "piglets", "lama glama", "ovis aries (domestic sheep)",
    "duck with tremor",
    "sus domesticus"]

livestock_isolation_sources = [
    "cattle", "chicken", "farm", "turkey feces", "wastewater from pig manure", "chicken feces",
    "pig", "swine nasal swab", "the intestine membrane of a diarrheic piglet",
    "intestine membrane of a diarrheic piglet", "pig manure",
    "isolated from the intestines of sick birds in the farm", "poultry establishment",
    "slaughtered pig", "poultry", "cat", "duck", "broiler", "pig faeces",
    "fecal samples from slaughtered sheep", "broiler chicken","chicken cecal content",
    "dog with mastitis", "goose anus swab",
    "porcine", "cattle hide", "goat; fecal","rumen fluid","pig environment","porcine feces",
    "cow raw milk", "cloacal swab",
    "lesion site (lung) of a dead turkey with colibacillosis", "laying hen withcolibacillosis",
    "goose faeces", "cattle slurry",
    "the respiratory tract of a pig with swine respiratory disease","cattle feces",
    "wool from pakistan", "chicken skin", "raw milk",
    "a nose swab sample of swine origin", "healthy broiler chicken","chicken manure","alpaca; fecal",
    "swine's gut", "swine final chilled carcass", "cow rumen", "cows", "dog",
    "bovine rumen", "cow dung","poultry carcass", "air of cow shed", "pig feces",
    "pig feed from feed plant", "fjerkrae",
    "tissue and/or biological fluid (swine)", "swine cecum", "raw cow milk","heifer vaginal mucus",
    "calf", "animal feed",
    "retail chicken","canine oral plaque","bovine pre-evisceration carcass at",
    "novine (bobby calf; fecal)", "a feces sample of chicken origin","lung of aborted horse fetus",
    "beef liver", "swine",
    "broiler chick cecum","isolated from the intestines of sick", "chicken cecum", "bulk pig ears",
    "cattle slurry", "cow feces",
    "ncsu equine educational unit", "cowshed of a farm", "a feces sample of swine origin",
    "sheep fecal sample", "digested slurry of dairy manure",
    "a anal swab sample of swine origin", "ear swab from dog", "goat", "healthy weaning piglets",
    "lesion site (lung) of a dead turkey", "chicken dung", "chicken trachea",
    "pooled sheep faecal samples collected from floor of farm", "poultry litter",
    "pooled pig faecal samples collected from floor of farm",
    "pooled cattle faecal samples collected from floor of farm",
    "rectal fecal grab samples from a commercial feedlot", "chicken tissues"]


anthropogenic_isolation_sources = [
    "cleaning system aquaculture", "clam larvae aquaculture", "mineral salt medium",
    "combined sewer", "oil", "roundup", "culture mutant", "turbot fish farm",
    "high concentration of fluoride", "fish larvae aquaculture", "sewage water",
    "combined sewer effluent", "feed additive imported from china", "wolfram mine tailing",
    "swab from a follow-up assessment of the tap handles and sink edges after a first disinfection attempt",
    "storage tank", "zc4 compost from a compost operation of sao paulo zoo",
    "swab from a hand-washing sink as part of the hospital routine surveillance program",
    "potash salt dump","bakery environment - hallway", "r2a medium", "polystyrene",
    "in the nicotine environment", "outdoor built environment", "farm kitchen",
    "clone isolated from the evolution experiment described in ketola et al. 2013 (doi: 10.1111/evo.12148)",
    "salt mine","bioreactor", "sewage & soil", "tung meal", "cooling tower",
    "surface of a polyethylene microplastic particle present in tank 6 of a marine aquarium containing stony-coral fragments and water maintained at 26 degree c", "aquarium water",
    "non-filtered water from the water column of tank 6 of a marine aquarium containing stony-coral fragments. water maintained at 26 degree c", "cheddar cheese factory", "mine",
    "telephone of nurse station","shower 3","geothermal power plant","sink aerator",
    "shelf","brewery environment", "lab strain", "coffee cup",
    "air in laboratory", "sp4 medium; universite de rennes", "industrial","knife", "wet market",
    "hospital", "open pond on an algae farm",
    "sf9 cell culture media", "laboratory","mold-colonized wall of an indoor",
    "power plant biotrickling filter",
    "shower 2", "room 7", "research laboratory","mycoplasma culture contaminant",
    "anaerobic digestion reactor",
    "heated-cooler unit water tank", "beer keg", "woodchip bioreactor",
    "wastewater treatment plant effluent",
    "bed sheets","automobile air-conditioning evaporator", "agricultural settling lagoon",
    "biocathode mcl",
    "polycyclic aromatic hydrocarbon", "agar plate", "brewery", "latex",
    "bakery environment - assembly production room",
    "anodic biofilm of glucose-fed microbial","milk powder production facility",
    "coke and gas plant treatment facilities",
    "rubber production plant territory", "laboratory stock","air conditioning system",
    "ensuite 7/8",
    "microbial mats from zloty stok gold mine","acid mine decant and tailings from uranium mine",
    "sedimentation pond in a zinc factory","acid mine drainage",
    "orfrc groundwater during biostimulation for uranium bioreduction",
    "sink drain", "chemical manufacturing sites","borehole hdn1; spa", "biofilm boat",
    "e-waste recycling site",
    "np30 nitrate/perchlorate-reducing mixed", "dust collector of pigpen",
    "acid mining effluent decantation pond",
    "anaerobic digester","cooling tower water", "dairy", "uninoculated hep-2 tissue culture",
    "tattoo ink",
    "arjo bathroom", "paper surface", "composted garbage","floor surface of biological laboratory",
    "hot water tap; geest office building", "industrial building air scrubber liquid",
    "iron water pipe", "sterile tools", "dental water pipeline","sungai pinang; penang; malaysia",
    "tar at the shipwreck site of tanker haven", "human septic tank","monitor panel",
    "western north pacific station s1", "ventilator", "coal-cleaning residues",
    "fresh-cut produce processing plant", "burns unit surveillance",
    "chilean kraft-pulp mill effluents", "automobile air conditional evaporator",
    "hospital sink", "oil-immersed sample from guaymas basin", "oil contaninated soil",
    "underground horizons of a flooded mine in russia",
    "red heat in salted hides (spoiled fish)", "colorado serum co.; anthrax spore vaccine",
    "electroporation of y pestis kim6+ in pcd1ap+", "solar panel array [envo:01000867]",
    "denitrifying; sulfide-oxidizing effluent-treatment plant",
    "high motility on swim plates; streptomycin resistant",
    "enrichment cultures from ucc lynda-stan",
    "base of single tree in sidewalk along rollins rd. evidence of many dogs; boston; mass",
    "monkey kidney tissue-culture fluids of the fh strain (eaton agent virus)",
    "conjugation assay", "lab", "culture maintained in leon; missing plasmid pscl3",
    "probitic products",
    "hygromycin b antibiotic bottle",
    "biofilm reactor", "nutrient broth"]

sediment_isolation_sources = [
    "mud","black mud", "hydrothermal vent area derived sediment", "envo:01001050",
    "marine mud", "muddy water", "pit mud", "tidal marsh", "freshwater mud",
    "mangrove swamp", "sea mud","mud flat", "the pit mud of a chinese liquor", "tidal flat sample",
    "muddy water",
    "mangrove wetland ecosystem", "subseafloor basaltic crust","lake mainaki salt",
    "mud from a salt lagoon"]

fungal_hosts = ["hymenochaete rubiginosa", "himantormia", "veronaeopsis simplex y34",
                "fusarium oxysporum f. sp. cucumerinum", "flammulina filiformis",
                "white-rot fungus phanerochaete chrysosporium"]

fungal_isolation_sources = ["mushroom substrate", "brewery yeast", "white-rot fungus phanerochaete",
                            "hyphae","white-rot fungus phanerochaete chrysosporium"]

human_isolation_sources = [
    "from patient with pneumonia", "blood", "stool","saliva","skin swab", "homo sapiens", "wound",
    "rectal swab", "skin", "fecal", "oral", "nosocomial environment",
    "human tissue biopsy","korean adult feces","groin","infant's throat",
    "sputum (cystic fibrosis patient)", "diarrheal patient",
    "clinical patient","intestine","respiratory tract","human clinical sample",
    "vomit of food poisoning patient", "faecal",
    "human skin swab","vomit from a food poisoning case",
    "china: the first affiliated hospital of zhejiang university","sputum",
    "human fecal sample","human blood","blood sample","urinary tract infection",
    "bronchial lavage","human feces", "anus swab",
    "wound infection","human","stool of a patient presenting with","breast abscess",
    "gastric","infant feces","eye",
    "human feces (woman; 24 years old)","human stool","urine","faeces","inflamed gingiva",
    "skin abscess",
    "tibia/osteomyelitis", "patient",
    "isolated from blood of a patient with","16-year-old boy with bubonic plague",
    "fecal sample","tissue",
    "blood of a hospitalized patient","brain","human feces (woman; 60 years old)",
    "foodborne disease surveillance",
    "korean infant feces", "healthy adult male gut","oropharynx","human foot","vein blood",
    "perirectal swab", "rectal", "anal swab",
    "blood specimen","rectal fecal grab samples from a","throat", "cystic fibrosis patient",
    "human skin tissue",
    "hospitalized patients","from patient with wound infection","human gi tract",
    "human case of meningitis", "supragingival plaque; periodontitis",
    "cancer patient stool", "milker's hand",
    "biopsy of antral stomach region from european patient with peptic ucler and chronic gastritis disease"]

na_hosts = ["environmental", "jamiecosley", "environment"]

na_isolation_sources = [
    "enviornmental","host's whole body","bacterial consortium", "missing", "unknown",
    "cell culture","afb-diseased colony","environment swab", "physical",
    "biological fluid","not available: to be reported later","environment",
    "r.j.roberts", "whole body", "atcc isolate", "feed",
    "isolate obtained from atcc:7955 nca","kimoto",
    "obscured","peter dedon mit","swab","dsm:3754","laboratory isolate","tsoundzou", "surface",
    "environmental surface", "ghana", "minnesota",
    "usa: wa: seattle","not available: to be reported later","not known","bone powder",
    "laboratory strain derived from cip","local geographical source", "dsmz isolate", "neb433", "bal",
    "air sacs","isolated from commercially available","contaminated wipes","protoplast breeding",
    "iam12617", "singapore", "ncimb", "valley",
    "swab with brown-gray powder",
    "anaerobic environments","environmental","neb269",
    "environmental swabs",
    "integrative microbiology research center; south china agricultural university (scau); tianhe district; guangzhou china",
    "coconi","csf","balf","china: zibo; shandong","wild type strain isolated from a natural source",
    "poland; warsaw area","agr",
    "Enviornmental"]

########################################################################################################################
#Script Begins Here
########################################################################################################################


## CRITICAL TODO: ~650 manual annotations are outdated, so many of them (~450) were marked NA but have been updated to have
## other values. Thus the code does not pick them up properly as it marks some cases as NA when they are not.
## for instance, some new annotations are annotated as "NA", probably due to an upstream bug in which the annotation
## fields changed in the original Genbank/GFF3 annotations. DOUBLE-CHECK ALL THESE CASES!

## For now, use manually-annotated-gbk-annotation-table.csv -- I stripped out the Manual_Annotation
## column to make a "pretend-input" table in order to test this code.

print("Annotation_Accession,host,isolation_source,Annotation")  ## print header

##with open("pretend-input-gbk-annotation-table.csv", "r") as annotation_fh: ##this version is just for debugging.
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
        elif isolation_source in na_isolation_sources and host == "na":
            annotation = "NA"
        elif host in na_hosts and isolation_source == "na":
            annotation = "NA"
        elif host in na_hosts and isolation_source in na_isolation_sources:
            annotation = "NA"

        ## Terrestrial
        elif host == "na" and isolation_source in terrestrial_isolation_sources:
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
        elif isolation_source in food_isolation_sources:
            annotation = "Food"
        elif (host == "na" and "milk" in isolation_source and "hand" not in isolation_source
              and "raw" not in isolation_source and "facility" not in isolation_source):
            annotation = "Food"
        elif host == "na" and "food" in isolation_source and "vomit" not in isolation_source and "sludge" not in isolation_source:
            annotation = "Food"
        elif host == "na" and "beef" in isolation_source:
            annotation = "Food"
        elif "kimchi" in isolation_source:
            annotation = "Food"
        elif host == "na" and "meat" in isolation_source and "facility" not in isolation_source:
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
        elif host == "sewage":
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
        elif host == "na" and " oil" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "facility" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "bioreactor" in isolation_source:
            annotation = "Anthropogenic-environment"
        elif host == "na" and "laboratory" in isolation_source:
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
        elif host == "na" and "human" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and "stool" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and "clinical" in isolation_source:
            annotation = "Human-host"
        elif host == "na" and isolation_source == "lung":
            annotation = "Human-host"

        print(','.join([annotation_accession, original_host, original_isolation_source, annotation]))

