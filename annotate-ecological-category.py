#!/usr/bin/env python

'''
annotate-ecological-category.py by Vincent Huang and Rohan Maddamsetti.

This script annotates the ecology for each analyzed genome, based on the metadata in the
'host' and 'isolation_source' tags in its genome file from Genbank.

Usage: python annotate-ecological-category.py > ../results/AR-gene-duplication/computationally-annotated-gbk-annotated-table.csv

'''

########################################################################################################################
## Keywords for annotating ecological categories.
########################################################################################################################
livestock_hosts = ["sus scrofa", "sus scrofa domesticus", "bos taurus", "piglet", "pig", "cattle","gallus gallus",
        "sheep placental tissue", "mouse-c57bl/6j", "swine", "chicken", "goose", "duck", "bovine","capra hircus",
        "gallus gallus domesticus", "porcine", "cow", "cat", "equus caballus ferus","pork", "sus","goat", "sheep",
        "equus ferus caballus", "bos bovis", "turkey", "rabbit","duck with tremor","bos primigenius taurus", "canine",
        "canis lupus familiaris", "bos taurus coreanae", "yak", "dog", "equine","ovis aries","young chicken",
        "dairy cow", "poultry", "piglets", "lama glama", "ovis aries (domestic sheep)"]

animal_hosts = ["donacia thalassina", "ixodes persulcatus", "tegillarca granosa", "edessa sp.","nematode",
        "tragelaphus strepsiceros", "apis mellifera", "botryllus sp.", "oncorhynchus kisutch","meleagris gallopavo",
        "neophocaena asiaeorientalis", "penaeus japonicus", "alligator", "macroplea mutica","tibetan antelope",
        "donacia provostii","donacia marginata", "ictalurus punctatus", "amblyomma cajennense","argas persicus","mink",
        "drosophila melanogaster oregon-r-modencode", "sablefish", "oncorhynchus mykiss", "atlantic salmon",
        "arvelius albopunctatus", "crassostrea gigas","phoca largha","amblyomma variegatum (cattle tick)","equus kiang",
        "penaeus vannamei","donacia versicolorea","nipponaphis monzeni","salmo salar","amazona sp.","plateau pika",
        "chroicocephalus novaehollandiae (australian silver","trachinotus ovatus","equus caballus", "shrimp", "crow",
        "galaxea fascicularis (stony coral)","bivalve mollusk","ondara zibethicus","aphis glycines","mus musculus",
        "donacia semicuprea","tuberolachnus salignus", "hippopotamus amphibius", "paralichthys olivaceus","bos mutus",
        "neoaliturus tenellus", "crucian carp", "ictalurus furcatus", "mercenaria mercenaria", "physopelta gutta",
        "apterostigma dentigerum","urocitellus undulatus", "mussels", "artemisaphis artemisicola", "aphis nerii","eel",
        "macroplea appendiculata","donacia clavipes","connochaetes taurinus","perinereis linea","trichonympha agilis",
        "drosophila melanogaster","dicentrarchus labrax (sea bass)", "rainbow trout", "donacia dentata", "honey bee",
        "ixodes scapularis","gut of wasp", "american cockroach", "trypoxylus dichotomus", "oyster larvae","parrot",
        "oreochromis niloticus", "ornithodoros hermsi", "oecophylla smaragdina", "euschistus servus","ixodes ricinus",
        "phakellia ventilabrum", "phoca vitulina", "megacopta punctatissima", "invertebrates","acyrthosiphon pisum",
        "donacia vulgaris","mouse", "mya arenaria oonogai makiyama", "nasonia vitripennis","cryptocercus punctulatus",
        "schizaphis graminum biotype i", "oreochromis", "cage-cultured red drum", "coral","ornithodoros sonrai", "fish",
        "gorilla", "hydra vulgaris aep", "amblyomma neumanni","plateumaris rustica","oryctes gigas","crab","honeybees",
        "syngnathus typhle","nezara viridula (cotton pathogen vector; southern", "litopenaeus vannamei","naegleria",
        "deinagkistrodon acutus", "aphis helianthi", "anguilla anguilla", "halyomorpha halys","antho dichotoma",
        "pthirus gorillae","marmota baibacina","ruddy shelduck", "marmota sibirica", "protaetia brevitarsis seulensis",
        "bactrocera oleae","rattus", "hyadaphis tataricae", "myzus persicae", "alces alces", "flatfish","toothfish",
        "acyrthosiphon kondoi", "sebastes schlegelii", "penaeus vannamei (shrimp)", "ovis","macrosiphum gaurae",
        "turbot", "pseudotrichonympha grassii (protist) in the gut of","cerambycidae sp.", "pediculus humanus corporis",
        "kangaroo", "cryptopsaras couesii","anas platyrhynchos", "penaeus setiferus", "apis mellifera mellifera",
        "camponotus chromaiodes", "caenorhabditis elegans","murgantia histrionica","ixodes pacificus","donacia proxima",
        "parakeet", "geronticus eremita","melanaphis sacchari","dendrolimus ibiricus","drosophila melanogaster oregon",
        "myzus persicae (green peach aphid)", "macrotermes barneyi", "glossina morsitans morsitans","nauphoeta cinerea",
        "oryctolagus cuniculus","corvus brachyrhynchos","drosophila melanogaster oregon-r modencode","apterostigma",
        "odocoileus virginianus","peromyscus leucopus", "cinara cedri", "brachycaudus cardui", "seabass", "rabbit",
        "plateumaris consimilis", "acanthamoeba polyphaga hn-3", "perameles bougainville","plateumaris braccata","rat",
        "hyalomma aegyptium","zebra", "apsterostigma", "baizongia pistaciae", "pediculus schaeffi", "sitobion avenae",
        "otospermophilus beecheyi (california ground","apodemus agrarius", "donacia fulgens", "bat", "skate",
        "paralichthys olivaceus (flounder)","dermacentor andersoni","donacia tomentosa","lacertilia","papio papio",
        "anomalops katoptron", "gull", "japanese eel", "salvelinus fontinalis","blattella germanica",
        "moschus berezovskii", "macaca silenus", "panulirus ornatus","pavona duerdeni", "odontobutis platycephala",
        "ostrea edulis (flat oyster)", "beluga whale", "shinkaia crosnieri","lutjanus guttatus (rose snapper)",
        "vultur gryphus","penaeus monodon","tanakia koreensis", "dicentrarchus labrax","phyllophaga sp.","bombyx mori",
        "zophobas atratus", "pigeon", "turtle", "penaeus (litopenaeus) vannamei (whiteleg shrimp)","donacia cinerea",
        "fulmars", "brevicoryne brassicae", "acanthamoeba", "blaberus giganteus","arctic char", "microlophium carnosum",
        "macaca mulatta","thelaxes californica","bactrocera dorsalis","apostichopus japonicus","ovis aries",
        "pthirus pubis","manis javanica (pangolin)", "sipalinus gigas", "geodia barretti","macrosiphoniella sanborni",
        "crassostrea gigas (pacific oyster)", "crassostrea virginica", "ellychnia corrusca","danaus plexippus",
        "glossina brevipalpis","marine sponge lissodendoryx isodictyalis in the","muscaphis stroyani","aphis urticata",
        "catfish", "morone chrysops x morone saxatilis", "pediculus humanus", "rattus norvegicus","plateumaris sericea",
        "bothriocroton concolor", "hyperomyzus lactucae","trimyema compressum", "macaca fascicularis", "diplonemea",
        "neohaemonia nigricornis", "haemaphysalis juxtakochi", "columba livia", "diaphorina citri", "prawn","deer",
        "pseudorca crassidens","acyrthosiphon lactucae","aphis nasturtii","giant panda", "cryptocercus kyebangensis",
        "japanese rhinoceros beetle larva","hydrophilus acuminatus", "discus", "hawk", "donacia sparganii", "oyster",
        "acinonyx jubatus","chaeturichthys stigmatias", "donacia cincticornis", "anas strepera","pan troglodytes",
        "aphis craccivora (cowpea aphid)", "calanoid copepod", "atlantic white-sided dolphin", "ctenocephalides felis",
        "pyropia yezoensis conchocelis", "pan paniscus", "lipaphis pseudobrassicae", "donacia piscatrix",
        "aphis craccivora","donacia crassipes", "donacia bicoloricornis", "donacia simplex", "draeculacephala minerva",
        "anodonta arcaeformis", "mastotermes darwiniensis", "toxic alexandrium minutum","coreoleuciscus splendidus",
        "platycercus elegans", "solea senegalensis","sanzinia madagascariensis volontany", "nezara viridula",
        "plateumaris pusilla","gilthead seabream", "ixodes pacificus (western blackleg tick)", "mus musculus c57bl/6j",
        "coregonus clupeaformis (lake whitefish)", "wax moth", "galleria mellonella","cyclopterus lumpus","cockatiel",
        "blatta orientalis", "seriola dumerili", "white stork", "canis latrans","drosophila neotestacea", "shrimps",
        "c57bl/6ntac mice", "phascolarctos cinereus", "neopsylla setosa", "haliclona simulans","ardea cinerea",
        "rhopalosiphum padi", "silurus asotus", "uroleucon ambrosiae","mustela putorius furo", "melanocetus johnsonii",
        "ailuropoda melanoleuca","spermophilus sp.", "squirrel"]

animal_isolation_sources = ["animal", "wild rat", "coho salmon", "fruit fly","diaphorina citri", "iguana",
        "efb-infected honey bee colony", "turban shell", "mussels", "rainbow trout","foulbrood of honeybees", "turbot",
        "marine sponge", "clam larvae aquaculture", "chamaeleonidae","ixodes spinipalpis", "barnacle at wood pile-on",
        "tilapia","fish larvae aquaculture", "wild pig; fecal", "shellfish","diarrheal snake diarrheal snake in",
        "insect larvae", "tissue; animal","yellowtail", "snake", "hydrophilus acuminatus", "crassostrea gigas","pigeon",
        "meleagris gallopavo","sea bass","mouse gut","avian","sockeye salmon","diseased labeo rohita fish","mealworm",
        "musca domestica","atlantic cod","siniperca scherzeri","midgut crypts of stink bug togo","konosirus punctatus",
        "crushed cell of callyspongia sp.", "tuna", "housefly", "bison", "coral","catfish", "mactra veneriformis",
        "shrimp", "black-headed gull","separated from the corpses of silkworms", "animal hide","tuna scrape; yellowfin",
        "aphis craccivora (cowpea aphid) on", "swiss alpine ibex feces","bearded dragon", "wild bird", "rattus rattus",
        "monkey kidney tissue-culture fluids of", "neophocaena phocaenoides", "loach","unidentified actinians",
        "biopsy of antral stomach region from", "avian","waxworms gut", "flea", "intestinal contents of termite"]

plant_hosts = ["actinidia","suaeda salsa","cladonia borealis","medicago truncatula",
        "allamanda cathartica","astragalus pelecinus","lebeckia ambigua","lotus corniculatus","trifolium repens",
        "datisca glomerata","beta vulgaris","robinia pseudoacacia","ficus religiosa l.","sporobolus anglicus",
        "gossypium hirsutum","leontopodium alpinum","tanacetum vulgare","catharanthus roseus","stereocaulon sp.",
        "agave americana l.","stachytarpheta glabra","elymus tsukushiensis","eucalyptus","androsace koso-poljanskii",
        "trifolium pratense","lemna minor","mimosa affinis","skeletonema marinoi strain st54","salsola stocksii",
        "biserrula pelecinus l.","kaempferia sp. rhizome","skeletonema marinoi strain ro5ac","masson pine","musa spp.",
        "nerium oleander","sugarcane","biserrula pelecinus","lichen","miscanthus giganteus","mimosa scabrella",
        "flower","mimosa flocculosa","trifolium uniflorum","eucalypti of eucalyptus","wild cotton",
        "trifolium spumosum l. (annual mediterranean","microcystis culture","medicago orbicularis","arabidopsis",
        "physcomitrella patens","cryptomeria fortunei","festuca arundinacea","oxytropis triphylla","usnea","rubus sp.",
        "achillea ptarmica l.","muscari","lespedeza cuneata","nerium oleander","tanacetum vulgare","eucommia ulmoides",
        "leiosporoceros dussii","oxytropis triphylla","vavilovia formosa","zantedeschia aethiopica","rosa",
        "leontopodium nivale","nusuttodinium aeruginosum","prosopis cineraria","oryza glumipatula","acacia farnesiana",
        "salvia splendens","panax ginseng", "vicia alpestris","phaeoceros","gossypium hirsutum",
        "commiphora wightii", "alhagi sparsifolia shap.", "acacia farnesiana","lobaria pulmonaria thallus", "lotus",
        "hyacinthus orientalis","sporobolus anglicus","malus prunifolia (crab apple)","tephrosia apollinea","musa sp.",
        "aeschynomene fluminensis","vavilovia formosa", "prunus sp.", "oenothera speciosa", "acaciella angustissima",
        "catharanthus roseus","dracaena sanderiana","parthenium argentatum gray (guayule shrubs)","medicago arborea",
        "euonymus sp.","eucalyptus grandis","broussonetia papyrifera","lichen stereocaulon sp.","lotus sp.","weed",
        "blasia pusilla","calliandra grandiflora","brittle root","digitaria eriantha","vavilovia formosa"]

plant_isolation_sources = ["flower of forsythia koreana", "plant rhizosphere", "peony rhizosphere",
        "algal phycosphere", "lotus corniculatus nodule","galega orientalis root nodule", "sugar beet rhizosphere",
        "flower of rhododendron schlippenbachii", "root of codonopsis pilosula",
        "chrysochromulina tobin phycosphere", "algae","artificially infected catharanthus",
        "peltigera membranacea thallus","angelica gigas nakai root surface",
        "plant","broussonetia papyrifera","lichen","seaweed; enteromorpha linza (coastal", "flower","leaf of fig tree",
        "plant roots","galega officinalis root nodule","halophyte rhizosphere","lichen thallus","ginseng","gladiolus",
        "sesbania spp. root nodule","plant leaf","origanum marjorana","malus sp.","algae","flower of shittah tree",
        "inner tissues of halophyte limonium","leaf","rhizosphere"]


agri_hosts = ["coffea", "grape", "citrus spp.", "prunus avium", "capsicum sp.","phaseolus vulgaris",
        "phaseolus", "punica granatum", "pyrus pyrifolia var. culta","alfalfa","cucumis sativus", "cucurbita maxima",
        "solanum melongena", "vitis vinifera","vitis vinifera l. cv. seto giants", "grapevine", "sakura tree", "rice",
        "chinese cabbage","kiwifruit", "pisum sativum l. (pea)", "common bean", "melon", "pepper", "glycine max",
        "actinidia deliciosa", "sugarcane", "actinidia chinensis", "lettuce", "olea europaea","grapefruit",
        "triticum aestivum", "areca catechu", "date palm", "zea mays l.", "ginger","plum", "cowpea", "potato", "wheat",
        "saccharum officinarum", "cotton", "soybean","mandarin orange", "fava bean", "capsicum annuum",
        "vitis vinifera l. cv. aurora black","ziziphus mauritiana lam", "pyrus pyrifolia", "strawberry",
        "sesame seedling", "tobacco","peanut", "carrot", "apple tree", "banana", "soybean (glycine max (l.) merrill)",
        "panicum miliaceum", "solanum lycopersicum", "plantain", "pear", "winter wheat","pyrus communis 'williams'",
        "actinida", "zea mays", "turfgrass", "prunus dulcis","glycine soja", "potato", "pepper plant", "sweet orange",
        "cicer arietinum","pisum sativum", "actinidia deliciosa 'hayward'", "arachis hypogaea","forage rape", "pear",
        "prunus cerasus (sour cherry)", "prunus dulcis", "sorghum bicolor", "phaseolus sp.","allium cepa",
        "vitis vinifera cv. 'izsaki sarfeher'","brassica rapa var. laciniifolia subvar.", "allium cepa l.",
        "gossypium hirsutum","solanum lycopersicoides", "vicia faba", "eggplant", "pear tree", "juglans regia",
                 "glycine max cv. jinju1", "lycopersicon esculentum", "camellia oleifera", "allium cepa (onion)"]

terrestrial_isolation_sources = ["iron hydroxide deposits", "solar salt", "air", "brazilian saline-akaline lake",
        "coal bed","moderate hot spring","salt marsh","alkaline pool submerged anode","kulunda steppe hypersaline lake",
        "solar saltern", "dust", "beach sand", "salar de atacama; atacama desert","weathered rock sample",
        "inside the caves of drach", "non-purified solar salt","gold-copper mine", "xinjiang aibi salt lake",
        "shar-burdin hypersaline soda lake", "marine solar saltern brine","soda lake magadi", "rock (ENVO:00001995)",
        "lava","cold-stored modified atmosphere", "brine", "saltern","colorado serum co.; anthrax spore",
        "taibei marine solar saltern near","orthoquartzite cave surface","cave","the bange salt-alkaline lake in tibet",
        "geothermal isolate", "salt crust", "coal spoil heap","flat; laminated microbial mat in a salt"]

soil_isolation_sources = ["agricultural field", "tundra wetland", "permafrost", "compost",
        "sphagnum peat from the bog obukhovskoe","paddy field; sungai manik; malaysia", "rice fields", "sphagnum peat",
        "sphagnum bog", "fertilizer", "ancient permafrost from mammoth","montane grasslands", "solar salt farm",
        "siberian permafrost","particulate matter", "pine forest", "mushroom compost", "manure compost"
        "composted cattle manure", "agricultural waste material","long-term organic manure fertilized"]

agri_isolation_sources = ["ucb-1 pistachio rootstock", "pineapple", "plants with moko disease",
                          "mangoes", "chives", "isolated from cucumber", "cabbage", "pea rhizosphere","almond drupe", "mexican lime leaf",
                          "healthy tomato plant","almond kernel (raw; variety carmel)", "the root of rice", "soybeans","wheat germ",
                          "sugarcane root", "originally isolated from olive trees","isolated from the tobacco substrate", "orange tree",
                          "hydroponic pots with potatoes", "ragi", "canola roots","leaf from rome apple cultivar", "cilantro","rice root",
                          "orange","maize stem", "strawberry leaf tissue", "rye silage", "fodder","stable grass silage","beans","carrot",
                          "corn","plant derived food stuff; onion; allium", "japanese pear","strawberry leaf tissue", "wheat root"]

agri_hosts = ["coffea", "grape", "citrus spp.", "prunus avium", "capsicum sp.","phaseolus vulgaris",
        "phaseolus", "punica granatum", "pyrus pyrifolia var. culta","alfalfa","cucumis sativus", "cucurbita maxima",
        "solanum melongena", "vitis vinifera","vitis vinifera l. cv. seto giants", "grapevine", "sakura tree", "rice",
        "chinese cabbage","kiwifruit", "pisum sativum l. (pea)", "common bean", "melon", "pepper", "glycine max",
        "actinidia deliciosa", "sugarcane", "actinidia chinensis", "lettuce", "olea europaea","grapefruit",
        "triticum aestivum", "areca catechu", "date palm", "zea mays l.", "ginger","plum", "cowpea", "potato", "wheat",
        "saccharum officinarum", "cotton", "soybean","mandarin orange", "fava bean", "capsicum annuum",
        "vitis vinifera l. cv. aurora black","ziziphus mauritiana lam", "pyrus pyrifolia", "strawberry",
        "sesame seedling", "tobacco","peanut", "carrot", "apple tree", "banana", "soybean (glycine max (l.) merrill)",
        "panicum miliaceum", "solanum lycopersicum", "plantain", "pear", "winter wheat","pyrus communis 'williams'",
        "actinida", "zea mays", "turfgrass", "prunus dulcis","glycine soja", "potato", "pepper plant", "sweet orange",
        "cicer arietinum","pisum sativum", "actinidia deliciosa 'hayward'", "arachis hypogaea","forage rape", "pear",
        "prunus cerasus (sour cherry)", "prunus dulcis", "sorghum bicolor", "phaseolus sp.","allium cepa",
        "vitis vinifera cv. 'izsaki sarfeher'","brassica rapa var. laciniifolia subvar.", "allium cepa l.",
        "gossypium hirsutum","solanum lycopersicoides", "vicia faba", "eggplant", "pear tree", "juglans regia",
        "glycine max cv. jinju1", "lycopersicon esculentum", "camellia oleifera", "allium cepa (onion)"]

marine_hosts = ["argopecten purpuratus", "seawater", "pyropia tenera", "trichodesmium erythraeum ims101",
        "rhodosorus marinus", "rhodobacter sphaeroides 2.4.1"]

marine_isolation_sources = ["salt lake", "hydrothermal vent area derived", "clam larvae aquaculture",
        "ocean", "isolated from the deepest ocean","chrysochromulina tobin phycosphere", "gulf of finland",
        "glass slide place on reef flat in","surface water of the southern north","antarctic iceberg",
        "hydrothermal vent", "acidic salty water","shallow tropical waters; normally from","open ocean water",
        "sulfidic waters (60 m) from the","shellfish hatchery", "bay of bengal", "pacific plankton",
        "zooplankton aquaculture", "culture of roseovarius sp. tm1035", "coastal","runway 10 reef (10-12 m)", "water",
        "microbial mat material from brackish","algae aquaculture", "hydrothermal precipitates",
        "shallow tropical waters; normally from", "beach sand","estuarine water [envo:01000301]", "150 km offshore",
        "cyanobacterial mat sample of antarctic","cyanobacterial aggregates","synechocystis sp. gt-l"]

freshwater_isolation_sources = ["creek","pond","rotifer aquaculture", "fish aquaculture","tianshan glacier",
        "fish larvae aquaculture", "glomma river","snow","microcystis culture","glacial stream", "algae aqualculture",
        "euglena gracilis from city ponds", "glacier from lahual spiti valley","river","pool at botanical garden"]

food_hosts = ["kimchi cabbage", "mung bean", "chicken", "kombucha scoby", "seujeot", "beef", "koumiss"]

food_isolation_sources = ["futsai", "koumiss", "soybean paste (chonggugjang)", "mashed potatoes",
        "light wheat beer", "sprouts", "3-day-old traditional semihard cheese","pork", "moto starter of sake",
        "degraded sugar thick juice","fresh alfalfa sprouts", "scallop", "chinese sausages", "ground turkey",
        "frozen peas", "palm wine", "natural whey culture from gruyere","homemade koumiss", "salad", "rotten eggs",
        "water kefir", "palm sap","peanut butter", "ground turkey", "sick cider", "onion","austrian hard cheese rind",
        "chinese pickle", "cheese", "raw mutton","nuruk; korean traditional beverage", "crab marinated in soy sauce",
        "honey", "ham", "bobby veal steak", "10' wieners", "natto", "yogurt","gochujang", "cheonggukjang",
        "traditional greek kasseri cheese","pork & cabbage dumplings", "lettuce", "new zealand cheese", "beer",
        "raw sausage", "beer contaminant", "scoby from kombucha tea", "calf liver","mustard pickles", "smoked fish",
        "doubanjiang", "retail turkey","fresh produce (lettuce; lollo bionda)", "blue berry", "apple cider","fried eel",
        "mac and cheese", "shelled pistachios", "marinated fish product","10 weeks old 45+ samso cheese", "salami",
        "mixed salads", "tibet kefir","sugar thick juice", "traditional dairy products", "sea bass", "ice cream",
        "water containing garlic and cabbage", "doenjang(soybean paste)","green chile ingredient",
        "koumiss (fermented mare's milk)", "pilsner beer","blown cheese", "pickled green chili peppers",
        "wheat sourdough", "cherry","pork barbeque", "vinegar", "dried tofu", "sauerkraut", "diced lamb", "salt",
        "korean traditional alcoholic beverage", "tibetan traditional dairy products","fish meal", "bagged lettuce",
        "apple juice from cider press", "sourdough","shrimp; frozen raw", "palm brown sugar",
        "tuna scrape; frozen yellowfin","isolate obtained from Campbell's Soup", "brine of stinky tofu", "mead",
        "french dry-type pork sausage", "sichuan pickle vegetables", "emmental","cheese production", "pork dumplings",
        "chinese egg powder", "stinky tofu","leaf vegetable", "retailed chichen", "home-made koumiss", "korean kefir"]

livestock_hosts = ["sus scrofa", "suxs scrofa domesticus", "bos taurus", "piglet", "pig", "cattle",
        "sheep placental tissue", "mouse-c57bl/6j", "swine", "chicken", "goose", "duck", "bovine","capra hircus",
        "gallus gallus", "gallus gallus domesticus", "porcine", "cow", "cat","equus caballus ferus", "pork", "sus",
        "goat", "sheep", "equus ferus caballus", "bos bovis","turkey", "rabbit", "bos primigenius taurus", "canine",
        "canis lupus familiaris","bos taurus coreanae", "yak", "dog", "equine", "ovis aries", "young chicken",
        "dairy cow", "poultry", "piglets", "lama glama", "ovis aries (domestic sheep)", "duck with tremor"]

livestock_isolation_sources = ["cattle", "meat processing facility", "chicken", "farm", "turkey feces",
        "feces", "fecal samples from slaughtered sheep", "broiler chicken","chicken cecal content", "dog with mastitis",
        "porcine", "skin","cattle hide", "goat; fecal","urine","rumen fluid","lung","pig environment","porcine feces",
        "lesion site(lung) of a dead turkey", "laying hen withcolibacillosis","goose faeces", "cattle slurry",
        "the respiratory tract of a pig with","cattle feces", "wool from pakistan", "chicken skin", "raw milk",
        "a nose swab sample of swine origin", "healthy broiler chicken","chicken manure","alpaca; fecal","swine's gut",
        "bovine rumen", "cow dung","poultry carcass", "air of cow shed", "pig feces", "pig feed from feed plant",
        "tissue and/or biological fluid (swine)", "swine cecum", "raw cow milk","heifer vaginal mucus", "calf",
        "red heat in salted hides (spoiled","retail chicken","canine oral plaque","bovine pre-evisceration carcass at",
        "novine (bobby calf; fecal)", "a feces sample of chicken origin","lung of aborted horse fetus", "beef liver",
        "broiler chick cecum","isolated from the intestines of sick", "chicken cecum", "bulk pig ears","cattle slurry"]

anthropogenic_isolation_sources = ["surface of a polyethylene microplastic", "cleaning system aquaculture",
        "mineral salt medium","milker's hand","clam larvae aquaculture","high concentration of fluoride",
        "potash salt dump","bakery environment - hallway", "r2a medium", "polystyrene", "salt mine","bioreactor",
        "telephone of nurse station","shower 3","geothermal power plant","sink aerator","shelf","brewery environment",
        "air in laboratory", "sp4 medium; universite de rennes", "industrial","knife", "wet market", "hospital",
        "sf9 cell culture media", "laboratory","mold-colonized wall of an indoor", "power plant biotrickling filter",
        "shower 2", "room 7", "research laboratory","mycoplasma culture contaminant", "anaerobic digestion reactor",
        "heated-cooler unit water tank", "beer keg", "woodchip bioreactor","wastewater treatment plant effluent",
        "bed sheets","automobile air-conditioning evaporator", "agricultural settling lagoon","biocathode mcl",
        "polycyclic aromatic hydrocarbon", "agar plate", "brewery", "latex","bakery environment - assembly",
        "anodic biofilm of glucose-fed microbial","milk powder production facility", "coke and gas plant treatment",
         "rubber production plant territory", "laboratory stock","air conditioning system", "ensuite 7/8",
        "microbial mats from zloty stok gold","acid mine decant and tailings from","hexachlorocyclohexane (hch)",
        "sedimentation pond in a zinc factory","acid mine drainage", "orfrc groundwater during biostimulation",
         "sink drain", "chemical manufacturing sites","borehole hdn1; spa", "biofilm boat", "e-waste recycling site",
         "np30 nitrate/perchlorate-reducing mixed", "dust collector of pigpen","acid mining effluent decantation pond",
        "anaerobic digester","cooling tower water", "dairy", "uninoculated hep-2 tissue culture", "tattoo ink",
        "arjo bathroom", "paper surface", "composted garbage","floor surface of biological laboratory",
        "meat processing facility","hot water tap; geest office building", "industrial building air scrubber",
        "iron water pipe", "sterile tools", "dental water pipeline","sungai pinang; penang; malaysia",
        "tar at the shipwreck site of tanker","undetoxicated lignocellulosic", "human septic tank","monitor panel",
        "western north pacific station s1", "ventilator", "coal-cleaning residues","fresh-cut produce processing plant",
        "chilean kraft-pulp mill effluents"]

sediment_isolation_sources = ["anaerobic top layer (5 to 10 cm) of", "mud","black mud",
        "microbial mat material from brackish", "marine mud", "muddy water", "pit mud", "tidal marsh", "freshwater mud",
        "mangrove swamp", "sea mud","mud flat", "the pit mud of a chinese liquor", "tidal flat sample","muddy water",
        "mangrove wetland ecosystem", "subseafloor basaltic crust","lake mainaki salt", "mud from a salt lagoon"]

fungal_hosts = ["hymenochaete rubiginosa", "himantormia", "veronaeopsis simplex y34","fusarium oxysporum f. sp. cucumerinum"]

fungal_isolation_sources = ["mushroom substrate", "brewery yeast", "white-rot fungus phanerochaete","hyphae"]

human_isolation_sources = ["from patient with pneumonia", "blood","feces","stool","saliva","skin swab",
        "human tissue biopsy","korean adult feces","groin","infant's throat","sputum (cystic fibrosis patient)",
        "clinical patient","intestine","respiratory tract","human clinical sample","vomit of food poisoning patient",
        "human skin swab","vomit from a food poisoning case","lung","china: the first affiliated hospital of","sputum",
        "human fecal sample","human blood","blood sample","urinary tract infection","bronchial lavage","human feces",
        "wound infection","human","stool of a patient presenting with","breast abscess","gastric","infant feces","eye",
        "human feces (woman; 24 years old)","human stool","urine","faeces","inflamed gingiva","skin abscess",
        "tibia/osteomyelitis","contaminated transfused red blood cell","human clinical sample (fatal pneumonic",
        "isolated from blood of a patient with","16-year-old boy with bubonic plague", "fecal sample","tissue",
        "blood of a hospitalized patient","brain","human feces (woman; 60 years old)","foodborne disease surveillance",
        "korean infant feces", "healthy adult male gut","oropharynx","human foot","vein blood", "perirectal swab",
        "blood specimen","rectal fecal grab samples from a","throat", "cystic fibrosis patient","human skin tissue",
        "hospitalized patients","from patient with wound infection","human gi tract","human case of meningitis"]

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

with open("../results/AR-gene-duplication/pretend-input-gbk-annotation-table.csv", "r") as annotation_fh: ## this version is just for debugging.
##with open("../results/AR-gene-duplication/gbk-annotation-table.csv", "r") as annotation_fh:
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
            
        ## Terrestrial
        elif host == "na" and isolation_source in terrestrial_isolation_sources:
            annotation = "Terrestrial"

        ## Soil
        elif host == "na" and "soil" in isolation_source and "contamina" not in isolation_source and "pollute" not in isolation_source:
            annotation = "Soil"
        elif host == "na" and isolation_source in soil_isolation_sources:
            annotation = "Soil"
        elif "soil" in host and "contamina" not in host and "pollute" not in host:
            annotation = "Soil"
        elif host == "na" and "rhizosphere" in isolation_source or "rhizosphere" in host:
            annotation = "Soil"

        ## Plant-host
        elif host not in agri_hosts and host in plant_hosts:
            annotation = "Plant-host"
        elif isolation_source in plant_isolation_sources:
            annotation = "Plant-host"

        ## Agriculture
        elif host not in plant_hosts and host in agri_hosts:
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

        ## Marine
        elif host in marine_hosts:
            annotation = "Marine"
        elif isolation_source in marine_isolation_sources:
            annotation = "Marine"
        elif "sea" in isolation_source:
            annotation = "Marine"
        elif "marine" in isolation_source:
            annotation = "Marine"

        ## Freshwater
        elif host == "rhodobacter sphaeroides 2.4.1":
            annotation = "Freshwater"
        elif "water" in isolation_source:
            annotation = "Freshwater"
        elif "spring" in isolation_source:
            annotation = "Freshwater"
        elif "lake" in isolation_source:
            annotation = "Freshwater"
        elif isolation_source in freshwater_isolation_sources:
            annotation = "Freshwater"
#
        ## Food
        elif host in food_hosts:
            annotation = "Food"
        elif host == "na" and isolation_source in food_isolation_sources:
            annotation = "Food"
        elif host == "na" and "milk" in isolation_source:
            annotation = "Food"
        elif host == "na" and "food" in isolation_source:
            annotation = "Food"
        elif host == "na" and "beef" in isolation_source:
            annotation = "Food"
        elif host == "na" and "kimchi" in isolation_source:
            annotation = "Food"
        elif host == "na" and "chicken" in isolation_source:
            annotation = "Food"
        elif host == "na" and "meat" in isolation_source:
            annotation = "Food"
        elif host == "na" and "ferment" in isolation_source:
            annotation = "Food"

        ## Livestock
        elif host in livestock_hosts:
            annotation = "Livestock"
        elif isolation_source in livestock_isolation_sources:
            annotation = "Livestock"

        ## Animal-host
        elif host in animal_hosts:
            annotation = "Animal-host"
        elif isolation_source in animal_isolation_sources:
            annotation = "Animal-host"

        ## Anthropogenic-environment
        elif host == "na" and isolation_source in anthropogenic_isolation_sources:
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

        ## Sediment
        elif host == "na" and "sediment" in isolation_source and "contamina" not in isolation_source and "pollute" not in isolation_source:
                annotation = "Sediment"
        elif host == "na" and isolation_source in sediment_isolation_sources:
            annotation = "Sediment"

        ## Fungal-host
        elif host in fungal_hosts or isolation_source in fungal_isolation_sources:
            annotation = "Fungal-host"

        ## Human-host
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and "homo" in host:
                annotation = "Human-host"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and host == "human":
            annotation = "Human-host"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and isolation_source in human_isolation_sources:
                annotation = "Human-host"

        print(','.join([annotation_accession, original_host, original_isolation_source, annotation]))
