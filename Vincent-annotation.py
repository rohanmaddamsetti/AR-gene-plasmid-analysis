#!/usr/bin/env python

## Annotation script by Vincent Huang.


print("Annotation_Accession,host,isolation_source,Manual_Annotation") ## print header
with open("../results/AR-gene-duplication/gbk-annotation-table.csv", "r") as annotation_fh:
    for line in annotation_fh:
        line = line.strip() ## remove trailing newline characters
        fields = line.split(',') ## split line into an array
        annotation_accession, host, isolation_source = fields
        manual_annotation = 'None'
        if host == "Homo sapiens":
            manual_annotation = "Human-host"
        elif host == "homo sapiens":
            manual_annotation = "Human-host"
        elif host == "Acyrthosiphon pisum":
            manual_annotation = "Animal-host"
        elif host == "Capsicum sp.":
            manual_annotation = "Agriculture"
        elif isolation_source == "olive fermentation":
            manual_annotation = "Agriculture"
        elif isolation_source == "dry beans":
            manual_annotation = "Agriculture"
        elif isolation_source == "production water from an oil well":
            manual_annotation = "Anthropogenic-environment"
        elif isolation_source == "clam larvae aquaculture":
            manual_annotation = "Marine"
        elif isolation_source == "forest soil":
            manual_annotation = "Soil"
        elif isolation_source == "air":
            manual_annotation = "Animal-host"
        elif isolation_source == "marine sediment":
            manual_annotation = "Sediment"
        elif host == "Citrus sinensis":
            manual_annotation = "Plant-host"
        elif isolation_source == "coal mine soil":
            manual_annotation = "Soil"
        elif isolation_source == "Baechu-kimchi":
            manual_annotation = "Plant-host"
        elif host == "alfalfa":
            manual_annotation = "Plant-host"
        elif isolation_source == "a nose swab sample of swine origin":
            manual_annotation = "Animal-host"
        elif isolation_source == "EFB-infected honey bee colony":
            manual_annotation = "Animal-host"
        elif isolation_source == "Buffalo milk":
            manual_annotation = "Animal-host"
        elif isolation_source == "floor surface of Biological laboratory":
            manual_annotation = "Anthropogenic-environment"
        elif isolation_source == "iron hydroxide deposits":
            manual_annotation = "Terrestrial"
        elif isolation_source == "fecal sample":
            manual_annotation = "Human-host"
        elif host == "Psoroma sp.":
            manual_annotation = "Plant-host"
        elif isolation_source == "hospitalized patients":
            manual_annotation = "Human-host"
        elif isolation_source == "sugarcane rhizosphere soil":
            manual_annotation = "Soil"
        elif host == "Apsterostigma":
            manual_annotation = "Animal-host"
        elif isolation_source == "sulfidic waters (60 m) from the Peruvian upwelling region":
            manual_annotation = "Marine"
        elif isolation_source == "Chinese sausages":
            manual_annotation = "Food"
        elif isolation_source == "Cowshed of a farm":
            manual_annotation = "Animal-host"
        elif host == "chicken":
            manual_annotation = "Animal-host"
        elif host == "pig":
            manual_annotation = "Animal-host"
        elif isolation_source == "sewage pipe":
            manual_annotation = "Anthropogenic-environment"
        elif isolation_source == "minced meat":
            manual_annotation = "Food"
        elif isolation_source == "poultry meat":
            manual_annotation = "Food"
        elif host == "chicken":
            manual_annotation = "Animal-host"
        elif host == "rice":
            manual_annotation = "Food"
        elif host == "Bird of prey":
            manual_annotation = "Animal-host"
        elif host == "Solanum lycopersicoides":
            manual_annotation = "Agriculture"
        elif host == "Salmo salar":
            manual_annotation = "Animal-host"
        elif host == "fish":
            manual_annotation = "Food"
        elif host == "pork":
            manual_annotation = "Food"
        elif isolation_source == "intestinal content":
            manual_annotation = "Animal-host"
        elif isolation_source == "stool of a patient presenting with scarlet-like fever in the Primorski region of the former USSR; in 1971 sent to the Institut Pasteur (Paris; France) by Dr. Timofeeva (Antiplague Institute; Irkoutsk)":
            manual_annotation = "Human-host"
        elif host == "Fusarium oxysporum f. sp. cucumerinum":
            manual_annotation = "Fungal-host"
        elif host == "Tibetan antelope":
            manual_annotation = "Animal-host"
        elif host == "Lotus corniculatus":
            manual_annotation = "Plant-host"
        elif host == "Musa sp.":
            manual_annotation = "Plant-host"
        elif host == "Bos primigenius taurus":
            manual_annotation = "Animal-host"
        elif host == "Penaeus vannamei":
            manual_annotation = "Animal-host"
        elif host == "Phaseolus":
            manual_annotation = "Soil"
        elif host == "Ixodes ricinus":
            manual_annotation = "Animal-host"
        elif host == "Donacia vulgaris":
            manual_annotation = "Animal-host"
        elif host == "Bos taurus":
            manual_annotation = "Animal-host"
        elif host == "honeybees":
            manual_annotation = "Animal-host"
        elif host == "dog":
            manual_annotation = "Animal-host"
        elif host == "eucalypti of Eucalyptus":
            manual_annotation = "Plant-host"
        elif host == "Gallus gallus":
            manual_annotation = "Animal-host"
        elif host == "Sus scrofa domesticus":
            manual_annotation = "Animal-host"
        elif host == "Papio papio":
            manual_annotation = "Animal-host"
        elif host == "citrus spp.":
            manual_annotation = "Plant-host"
        elif host == "Macroplea appendiculata":
            manual_annotation = "Animal-host"
        elif host == "Solanum tuberosum":
            manual_annotation = "Plant-host"
        elif host == "Sus scrofa":
            manual_annotation = "Animal-host"
        elif host == "Catharanthus roseus":
            manual_annotation = "Plant-host"
        elif host == "Equus caballus":
            manual_annotation = "Animal-host"
        elif host == "Pepper":
            manual_annotation = "Plant-host"
        elif host == "sugarcane":
            manual_annotation = "Plant-host"
        elif host == "Drosophila melanogaster":
            manual_annotation = "Animal-host"
        elif host == "Protaetia brevitarsis seulensis":
            manual_annotation = "Animal-host"
        elif host == "Apis mellifera":
            manual_annotation = "Animal-host"
        elif host == "Gossypium sp.":
            manual_annotation = "Plant-host"
        elif host == "Medicago sativa L. subsp. ambigua":
            manual_annotation = "Plant-host"
        elif host == "common bean":
            manual_annotation = "Plant-host"
        elif host == "Syngnathus typhle":
            manual_annotation = "Terrain"
        elif host == "Wheat":
            manual_annotation = "Plant-host"
        elif host == "red alga":
            manual_annotation = "Marine"
        elif host == "Goat":
            manual_annotation = "Animal-host"
        elif host == "Penaeus japonicus":
            manual_annotation = "Animal-host"
        elif isolation_source == "leaf":
            manual_annotation = "Plant-host"
        elif host == "Beef":
            manual_annotation = "Animal-host"
        elif host == "flower":
            manual_annotation = "Plant-host"
        elif host == "Wild boar":
            manual_annotation = "Animal-host"
        elif host == "Atlantic white-sided dolphin":
            manual_annotation = "Animal-host"
        elif host == "Geodia barretti":
            manual_annotation = "Animal-host"
        elif host == "Canis lupus familiaris":
            manual_annotation = "Animal-host"
        elif host == "Elymus tsukushiensis":
            manual_annotation = "Soil"
        elif host == "porcine":
            manual_annotation = "Animal-host"
        elif host == "Soybean (Glycine max (L.) Merrill)":
            manual_annotation = "Plant-host"
        elif host == "Turbot":
            manual_annotation = "Animal-host"
        elif host == "Androsace koso-poljanskii":
            manual_annotation = "Plant-host"
        elif isolation_source == "flower":
            manual_annotation = "Plant-host"
        elif host == "duck":
            manual_annotation = "Animal-host"
        elif host == "Triticum aestivum L.":
            manual_annotation = "Plant-host"
        elif host == "Halyomorpha halys":
            manual_annotation = "Animal-host"
        elif host == "mouse":
            manual_annotation = "Animal-host"
        elif host == "Citrus paradisi":
            manual_annotation = "Plant-host"
        elif host == "Columba livia":
            manual_annotation = "Animal-host"
        elif host == "Sesame":
            manual_annotation = "Plant-host"
        elif host == "banana":
            manual_annotation = "Plant-host"
        elif host == "Oryza glumipatula":
            manual_annotation = "Plant-host"
        elif host == "porcine":
            manual_annotation = "Animal-host"
        elif isolation_source == "seed":
            manual_annotation = "Plant-host"
        elif host == "Neophocaena asiaeorientalis":
            manual_annotation = "Animal-host"
        elif host == "Citrus hystrix":
            manual_annotation = "Plant-host"
        elif host == "Canis lupus familiaris":
            manual_annotation = "Animal-host"
        elif isolation_source == "Beetle Malpighian tubules":
            manual_annotation = "Animal-host"
        elif host == "Aphis urticata":
            manual_annotation = "Animal-host"
        elif host == "poultry":
            manual_annotation = "Animal-host"
        elif isolation_source == "larva":
            manual_annotation = "Animal-host"
        elif host == "Broussonetia papyrifera":
            manual_annotation = "Plant-host"
        elif host == "Bos grunniens":
            manual_annotation = "Animal-host"
        elif host == "Bos taurus coreanae":
            manual_annotation = "Animal-host"
        elif host == "Skeletonema marinoi strain RO5AC":
            manual_annotation = "Sediment"
        elif host == "Diplonemea":
            manual_annotation = "Animal-host"
        elif host == "Macrotermes barneyi":
            manual_annotation = "Animal-host"
        elif host == "Glycine max":
            manual_annotation = "Plant-host"
        elif host == "Macrotermes barneyi":
            manual_annotation = "Animal-host"
        elif host == "Arabidopsis":
            manual_annotation = "Plant-host"
        elif host == "Apple":
            manual_annotation = "Plant-host"
        elif host == "Pyrus communis 'Williams'":
            manual_annotation = "Plant-host"
        elif host == "Glossina morsitans morsitans":
            manual_annotation = "Plant-host"
        elif host == "Bovine":
            manual_annotation = "Plant-host"
        elif isolation_source == "kidney":
            manual_annotation = "Animal-host"
        elif host == "Skeletonema marinoi strain RO5AC":
            manual_annotation = "Sediment"
        elif host == "Cotinus coggygria":
            manual_annotation = "Plant-host"
        elif host == "coral":
            manual_annotation = "Marine"
        elif host == "Sanzinia madagascariensis volontany":
            manual_annotation = "Animal-host"
        elif host == "Achillea ptarmica L.":
            manual_annotation = "Plant-host"
        elif host == "Physopelta gutta":
            manual_annotation = "Animal-host"
        elif host == "Oncorhynchus mykiss":
            manual_annotation = "Animal-host"
        elif isolation_source == "rhizosphere":
            manual_annotation = "Soil"
        elif host == "Trifolium pratense":
            manual_annotation = "Plant-host"
        elif isolation_source == "choana":
            manual_annotation = "Animal-host"
        elif host == "Giant panda":
            manual_annotation = "Animal-host"
        elif host == "cat":
            manual_annotation = "Animal-host"
        elif host == "Tanacetum vulgare":
            manual_annotation = "Plant-host"
        elif host == "mouse-C57Bl/6J":
            manual_annotation = "Animal-host"
        elif host == "Ondatra zibethicus":
            manual_annotation = "Animal-host"
        elif host == "Sipalinus gigas":
            manual_annotation = "Animal-host"
        elif host == "Crassostrea gigas (Pacific oyster)":
            manual_annotation = "Animal-host"
        elif isolation_source == "milk":
            manual_annotation = "Animal-host"
        elif host == "Oncorhynchus kisutch":
            manual_annotation = "Marine"
        elif host == "kiwifruit":
            manual_annotation = "Plant-host"
        elif host == "Pthirus gorillae":
            manual_annotation = "Animal-host"
        elif host == "Mya arenaria oonogai Makiyama":
            manual_annotation = "Animal-host"
        elif host == "Glycine max cv. Jinju1":
            manual_annotation = "Plant-host"
        elif host == "Argas persicus":
            manual_annotation = "Animal-host"
        elif host == "Alvinocaris longirostris":
            manual_annotation = "Animal-host"
        elif host == "Capra hircus":
            manual_annotation = "Animal-host"
        elif host == "Tephrosia apollinea":
            manual_annotation = "Plant-host"
        elif host == "Trypoxylus dichotomus":
            manual_annotation = "Animal-host"
        elif isolation_source == "rhizosphere":
            manual_annotation = "Soil"
        elif host == "Trachinotus ovatus":
            manual_annotation = "Animal-host"
        elif host == "Canis familiaris":
            manual_annotation = "Animal-host"
        elif isolation_source == "spleen":
            manual_annotation = "Animal-host"
        elif host == "Oryctolagus cuniculus":
            manual_annotation = "Animal-host"
        elif host == "Galleria mellonella":
            manual_annotation = "Animal-host"
        elif isolation_source == "plant":
            manual_annotation = "Plant-host"
        elif host == "Phaseolus vulgaris":
            manual_annotation = "Plant-host"
        elif host == "Oncorhynchus kisutch":
            manual_annotation = "Animal-host"
        elif host == "Lettuce":
            manual_annotation = "Plant-host"
        elif host == "Bivalve mollusk":
            manual_annotation = "Marine"
        elif host == "Oncorhynchus mykiss":
            manual_annotation = "Animal-host"
        elif host == "Solea senegalensis":
            manual_annotation = "Animal-host"
        elif host == "wax moth":
            manual_annotation = "Animal-host"
        elif isolation_source == "sponge":
            manual_annotation = "Animal-host"
        elif host == "mung bean":
            manual_annotation = "Plant-host"
        elif host == "Homo sapiens newborn":
            manual_annotation = "Human-host"
        elif host == "Anas platyrhynchos":
            manual_annotation = "Animal-host"
        elif host == "Aphis craccivora (cowpea aphid)":
            manual_annotation = "Animal-host"
        elif host == "Schizaphis graminum biotype I":
            manual_annotation = "Animal-host"
        elif host == "Penaeus vannamei":
            manual_annotation = "Animal-host"
        elif host == "Mustela putorius furo":
            manual_annotation = "Animal-host"
        elif host == "Leontopodium alpinum":
            manual_annotation = "Terrain"
        elif host == "Pyropia yezoensis conchocelis":
            manual_annotation = "Marine"
        elif host == "Miscanthus giganteus":
            manual_annotation = "Plant-host"
        elif host == "Prunus avium":
            manual_annotation = "Plant-host"
        elif host == "crucian carp":
            manual_annotation = "Animal-host"
        elif isolation_source == "lung":
            manual_annotation = "Animal-host"
        elif isolation_source == "water":
            manual_annotation = "Animal-host"
        elif isolation_source == "nodule":
            manual_annotation = "Plant-host"
        elif host == "egret":
            manual_annotation = "Animal-host"
        elif host == "Pyropia":
            manual_annotation = "Plant-host"
        elif host == "Rainbow trout":
            manual_annotation = "Animal-host"
        elif host == "Moschus berezovskii":
            manual_annotation = "Animal-host"
        elif host == "Meleagris gallopavo":
            manual_annotation = "Animal-host"
        elif host == "Hyalomma aegyptium":
            manual_annotation = "Animal-host"
        elif host == "Mus musculus":
            manual_annotation = "Animal-host"
        elif host == "Androsace koso-poljanskii":
            manual_annotation = "Plant-host"
        elif host == "Phaseolus vulgaris":
            manual_annotation = "Plant-host"
        elif host == "Salsola stocksii":
            manual_annotation = "Plant-host"
        elif host == "weed":
            manual_annotation = "Plant-host"
        elif host == "flatfish":
            manual_annotation = "Animal-host"
        elif host == "Ovis aries":
            manual_annotation = "Animal-host"
        elif host == "Argopecten purpuratus":
            manual_annotation = "Animal-host"
        elif host == "Homo sapiens (female)":
            manual_annotation = "Human-host"
        elif host == "Gossypium hirsutum":
            manual_annotation = "Agriculture"
        elif host == "Citrus aurantifolia":
            manual_annotation = "Plant-host"
        elif host == "Amazona sp.":
            manual_annotation = "Animal-host"
        elif host == "Sporobolus anglicus":
            manual_annotation = "Animal-host"
        elif host == "Calliandra grandiflora":
            manual_annotation = "Plant-host"
        elif host == "Haemaphysalis juxtakochi":
            manual_annotation = "Plant-host"
        elif host == "Pisum sativum":
            manual_annotation = "Plant-host"
        elif host == "Equus ferus caballus":
            manual_annotation = "Animal-host"
        elif host == "Cyclopterus lumpus":
            manual_annotation = "Animal-host"
        elif host == "Avena sativa":
            manual_annotation = "Plant-host"
        elif host == "Cryptopsaras couesii":
            manual_annotation = "Animal-host"
        elif host == "Allium cepa":
            manual_annotation = "Plant-host"
        elif host == "Ornithodoros hermsi":
            manual_annotation = "Animal-host"
        elif host == "Brassica rapa subsp. pekinensis":
            manual_annotation = "Plant-host"
        elif host == "Peromyscus leucopus":
            manual_annotation = "Animal-host"
        elif host == "Citrus aurantifolia":
            manual_annotation = "Plant-host"
        elif host == "Microlophium carnosum":
            manual_annotation = "Animal-host"
        elif host == "sludge":
            manual_annotation = "Anthropogenic-environment"
        elif host == "crow":
            manual_annotation = "Animal-host"
        elif host == "winter wheat":
            manual_annotation = "Plant-host"
        elif host == "Cladonia borealis":
            manual_annotation = "Plant-host"
        elif host == "Olea europaea":
            manual_annotation = "Plant-host"
        elif host == "invertebrates":
            manual_annotation = "Soil"
        elif host == "Oncorhynchus mykiss":
            manual_annotation = "Animal-host"
        elif host == "Glycine soja":
            manual_annotation = "Plant-host"
        elif host == "Penaeus japonicus":
            manual_annotation = "Animal-host"
        elif host == "Phoca vitulina":
            manual_annotation = "Animal-host"
        elif host == "Bos taurus coreanae":
            manual_annotation = "Animal-host"
        elif host == "Phaseolus":
            manual_annotation = "Soil"
        elif host == "Glycine max L. Merr.":
            manual_annotation = "Plant-host"
        elif host == "Bos primigenius taurus":
            manual_annotation = "Animal-host"
        elif host == "pheasant":
            manual_annotation = "Animal-host"
        elif host == "Oryza sativa":
            manual_annotation = "Plant-host"
        elif host == "Crassostrea gigas":
            manual_annotation = "Animal-host"
        elif host == "Solanum melongena":
            manual_annotation = "Plant-host"
        elif host == "Allamanda cathartica":
            manual_annotation = "Plant-host"
        elif host == "Melanaphis sacchari":
            manual_annotation = "Animal-host"
        elif host == "gram-positive bacteria":
            manual_annotation = "Animal-host"
        elif host == "Mus musculus C57Bl/6J":
            manual_annotation = "Animal-host"
        elif host == "Panax ginseng":
            manual_annotation = "Plant-host"
        elif host == "Bothriocroton concolor":
            manual_annotation = "Animal-host"
        elif host == "Allium cepa L.":
            manual_annotation = "Plant-host"
        elif host == "common bean":
            manual_annotation = "Plant-host"
        elif host == "Parthenium argentatum Gray (guayule shrubs)":
            manual_annotation = "Animal-host"
        elif host == "fulmars":
            manual_annotation = "Animal-host"
        elif host == "Marmota himalayana":
            manual_annotation = "Animal-host"
        elif host == "Citrus limon":
            manual_annotation = "Plant-host"
        elif host == "Crassostrea gigas":
            manual_annotation = "Animal-host"
        elif host == "Neophocaena asiaeorientalis":
            manual_annotation = "Animal-host"
        elif host == "Pisum sativum l. (pea)":
            manual_annotation = "Plant-host"
        print(','.join([annotation_accession, host, isolation_source, manual_annotation]))
