### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 86e1ee96-f314-11ea-03f6-0f549b79e7c9
begin
	using Pkg
	Pkg.activate(mktempdir())
end

# ╔═╡ a4937996-f314-11ea-2ff9-615c888afaa8
begin
	Pkg.add([
			"Compose",
			"Colors",
			"PlutoUI",
			"DataFrames",
			"DataFramesMeta",
			"CSVFiles",
			])

	using Colors
	using PlutoUI
	using Compose
	using LinearAlgebra
	using DataFrames
	using DataFramesMeta
	using CSVFiles
end

# ╔═╡ 85cfbd10-f384-11ea-31dc-b5693630a4c5
md"""
# **ngrams.jl**: _N-gram analysis of duplicated genes in bacterial genomes_

by Rohan Maddamsetti. 

Cluster functional annotation from duplicated genes associated with different environmental metadata.  


"""

# ╔═╡ 938185ec-f384-11ea-21dc-b56b7469f798
md"_Let's create a package environment:_"

# ╔═╡ 3cc688d2-f996-11ea-2a6f-0b4c7a5b74c2
md"""
Let's make a **frequency table** for n-grams in bacterial genomes.  

Let's import the duplicated gene data, the annotations, and then join the annotations to the duplicated gene data.

$(html"<br>")

Remember that the $(html"<img src='https://cdn.jsdelivr.net/gh/ionic-team/ionicons@5.0.0/src/svg/eye-outline.svg' style='width: 1em; height: 1em; margin-bottom: -.2em;'>") button on the left of a cell will show or hide the code.

"""

# ╔═╡ 829eea8e-05f0-11eb-2a9f-cf1b0c339053
begin
	duplicated_genes_df = DataFrame(load("../results/AR-gene-duplication/duplicate-proteins.csv"))
	annotation_df = DataFrame(load("../data/manually-curated-gbk-annotation-table.csv"))
	annotated_duplicates_df = leftjoin(duplicated_genes_df,annotation_df,on=:Annotation_Accession)
end

# ╔═╡ 28195e06-05f4-11eb-3997-2377c3e3b124
names(annotated_duplicates_df)

# ╔═╡ 9424c438-0652-11eb-2329-4b90f2f7cecf
md""" Filter all duplicated proteins that match MGE keywords.
"""

# ╔═╡ ae61d73c-0652-11eb-38b5-df004da5eb69
IS_keywords = r"IS|transposon|Transposase|transposase|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin"

# ╔═╡ b6e8e120-0652-11eb-3af6-dbc37e3ae9b5
filtered_annotated_duplicates_df = @where(annotated_duplicates_df, occursin(IS_keywords, :product))

# ╔═╡ 0d3aaa6a-0651-11eb-0034-d77db75bb34c
annotated_duplicates_df[:product]

# ╔═╡ a094e2ac-f92d-11ea-141a-3566552dd839
md"""
#### _Data cleaning_
Looking at the sample, we see that it is quite _messy_: it may contain punctuation, accented letters and numbers. For our analysis, we are only interested in lower-case words. We are going to clean the data using the Julia function `filter`. 
"""

# ╔═╡ faf80910-0650-11eb-0798-dd4b1628700e
alphabet = ['a':'z' ; ' ']   # includes the space character

# ╔═╡ 5c74a052-f92e-11ea-2c5b-0f1a3a14e313
function isinalphabet(character)
	character ∈ alphabet
end

# ╔═╡ 05f0182c-f999-11ea-0a52-3d46c65a049e
md"""
$(html"<br>")

We are not interested in the case of letters (i.e. `'A'` vs `'a'`), so we want to map these to lower case _before_ we apply our filter. If we don't, all upper case letters would get deleted.
"""

# ╔═╡ 24860970-fc48-11ea-0009-cddee695772c
import Unicode

# ╔═╡ 734851c6-f92d-11ea-130d-bf2a69e89255
"""
Turn `"áéíóúüñ asdf"` into `"aeiouun asdf"`.
"""
unaccent(str) = Unicode.normalize(str, stripmark=true)

# ╔═╡ d67034d0-f92d-11ea-31c2-f7a38ebb412f
samples = (
	English = """
Although the word forest is commonly used, there is no universally recognised precise definition, with more than 800 definitions of forest used around the world.[4] Although a forest is usually defined by the presence of trees, under many definitions an area completely lacking trees may still be considered a forest if it grew trees in the past, will grow trees in the future,[9] or was legally designated as a forest regardless of vegetation type.[10][11]
	
The word forest derives from the Old French forest (also forès), denoting "forest, vast expanse covered by trees"; forest was first introduced into English as the word denoting wild land set aside for hunting[14] without the necessity in definition of having trees on the land.[15] Possibly a borrowing, probably via Frankish or Old High German, of the Medieval Latin foresta, denoting "open wood", Carolingian scribes first used foresta in the Capitularies of Charlemagne specifically to denote the royal hunting grounds of the King. The word was not endemic to Romance languages, e. g. native words for forest in the Romance languages derived from the Latin silva, which denoted "forest" and "wood(land)" (confer the English sylva and sylvan); confer the Italian, Spanish, and Portuguese selva; the Romanian silvă; and the Old French selve, and cognates in Romance languages, e. g. the Italian foresta, Spanish and Portuguese floresta, etc., are all ultimately derivations of the French word. 
""",
	Spanish =  """
Un bosque es un ecosistema donde la vegetación predominante la constituyen los árboles y matas.1​ Estas comunidades de plantas cubren grandes áreas del globo terráqueo y funcionan como hábitats para los animales, moduladores de flujos hidrológicos y conservadores del suelo, constituyendo uno de los aspectos más importantes de la biosfera de la Tierra. Aunque a menudo se han considerado como consumidores de dióxido de carbono atmosférico, los bosques maduros son prácticamente neutros en cuanto al carbono, y son solamente los alterados y los jóvenes los que actúan como dichos consumidores.2​3​ De cualquier manera, los bosques maduros juegan un importante papel en el ciclo global del carbono, como reservorios estables de carbono y su eliminación conlleva un incremento de los niveles de dióxido de carbono atmosférico.

Los bosques pueden hallarse en todas las regiones capaces de mantener el crecimiento de árboles, hasta la línea de árboles, excepto donde la frecuencia de fuego natural es demasiado alta, o donde el ambiente ha sido perjudicado por procesos naturales o por actividades humanas. Los bosques a veces contienen muchas especies de árboles dentro de una pequeña área (como la selva lluviosa tropical y el bosque templado caducifolio), o relativamente pocas especies en áreas grandes (por ejemplo, la taiga y bosques áridos montañosos de coníferas). Los bosques son a menudo hogar de muchos animales y especies de plantas, y la biomasa por área de unidad es alta comparada a otras comunidades de vegetación. La mayor parte de esta biomasa se halla en el subsuelo en los sistemas de raíces y como detritos de plantas parcialmente descompuestos. El componente leñoso de un bosque contiene lignina, cuya descomposición es relativamente lenta comparado con otros materiales orgánicos como la celulosa y otros carbohidratos. Los bosques son áreas naturales y silvestre 
""" |> unaccent,
)

# ╔═╡ 8d3bc9ea-f9a1-11ea-1508-8da4b7674629
md"""
$(html"<br>")

👉 Let's put everything together. Write a function `clean` that takes a string, and returns a _cleaned_ version, where:
- accented letters are replaced by their base characters;
- upper-case letters are converted to lower case;
- it is filtered to only contain characters from `alphabet`
"""

# ╔═╡ 4affa858-f92e-11ea-3ece-258897c37e51
function clean(text)
	# we turn everything to lowercase to keep the number of letters small
	return filter(isinalphabet, lowercase(unaccent(text)))
end

# ╔═╡ 2680b506-f9a3-11ea-0849-3989de27dd9f
first_sample = clean(first(samples))

# ╔═╡ 571d28d6-f960-11ea-1b2e-d5977ecbbb11
function letter_frequencies(txt)
	ismissing(txt) && return missing
	f = count.(string.(alphabet), txt)
	f ./ sum(f)
end

# ╔═╡ 6a64ab12-f960-11ea-0d92-5b88943cdb1a
sample_freqs = letter_frequencies(first_sample)

# ╔═╡ 603741c2-f9a4-11ea-37ce-1b36ecc83f45
md"""
The result is a 27-element array, with values between `0.0` and `1.0`. These values correspond to the _frequency_ of each letter. 

`sample_freqs[i] == 0.0` means that the $i$th letter did not occur in your sample, and 
`sample_freqs[i] == 0.1` means that 10% of the letters in the sample are the $i$th letter.

To make it easier to convert between a character from the alphabet and its index, we have the following function:
"""

# ╔═╡ b3de6260-f9a4-11ea-1bae-9153a92c3fe5
index_of_letter(letter) = findfirst(isequal(letter), alphabet)

# ╔═╡ a6c36bd6-f9a4-11ea-1aba-f75cecc90320
index_of_letter('a'), index_of_letter('b'), index_of_letter(' ')

# ╔═╡ 6d3f9dae-f9a5-11ea-3228-d147435e266d
md"""
$(html"<br>")

👉 Which letters from the alphabet did not occur in the sample?
"""

# ╔═╡ 92bf9fd2-f9a5-11ea-25c7-5966e44db6c6
unused_letters = let
	['j','q','z']
end

# ╔═╡ 01215e9a-f9a9-11ea-363b-67392741c8d4
md"""
**Random letters at the correct frequencies:**
"""

# ╔═╡ 8ae13cf0-f9a8-11ea-3919-a735c4ed9e7f
md"""
By considering the _frequencies_ of letters in English, we see that our model is already a lot better!

Our next observation is that **some letter _combinations_ are more common than others**. Our current model thinks that `potato` is just as 'English' as `ooaptt`. In the next section, we will quantify these _transition frequencies_, and use it to improve our model.
"""

# ╔═╡ 343d63c2-fb58-11ea-0cce-efe1afe070c2


# ╔═╡ b5b8dd18-f938-11ea-157b-53b145357fd1
function rand_sample(frequencies)
	x = rand()
	findfirst(z -> z >= x, cumsum(frequencies ./ sum(frequencies)))
end

# ╔═╡ 0e872a6c-f937-11ea-125e-37958713a495
function rand_sample_letter(frequencies)
	alphabet[rand_sample(frequencies)]
end

# ╔═╡ fbb7c04e-f92d-11ea-0b81-0be20da242c8
function transition_counts(cleaned_sample)
	[count(string(a, b), cleaned_sample)
		for a in alphabet,
			b in alphabet]
end

# ╔═╡ 80118bf8-f931-11ea-34f3-b7828113ffd8
normalize_array(x) = x ./ sum(x)

# ╔═╡ 7f4f6ce4-f931-11ea-15a4-b3bec6a7e8b6
transition_frequencies = normalize_array ∘ transition_counts;

# ╔═╡ d40034f6-f9ab-11ea-3f65-7ffd1256ae9d
transition_frequencies(first_sample)

# ╔═╡ 689ed82a-f9ae-11ea-159c-331ff6660a75
md"What we get is a **27 by 27 matrix**. Each entry corresponds to a character pair. The _row_ corresponds to the first character, the _column_ is the second character. Let's visualize this:"

# ╔═╡ 0b67789c-f931-11ea-113c-35e5edafcbbf
md"""
Answer the following questions with respect to the **cleaned English sample text**, which we called `first_sample`. Let's also give the transition matrix a name:
"""

# ╔═╡ 6896fef8-f9af-11ea-0065-816a70ba9670
sample_freq_matrix = transition_frequencies(first_sample);

# ╔═╡ 39152104-fc49-11ea-04dd-bb34e3600f2f
if first_sample === missing
	md"""
	!!! danger "Don't worry!"
	    👆 These errors will disappear automatically once you have completed the earlier exercises!
	"""
end

# ╔═╡ 2f8dedfc-fb98-11ea-23d7-2159bdb6a299
md"""
We can use the measured transition frequencies to generate text in a way that it has **the same transition frequencies** as our original sample. Our generated text is starting to look like real language!
"""

# ╔═╡ b7446f34-f9b1-11ea-0f39-a3c17ba740e5
@bind ex23_sample Select([v => String(k) for (k,v) in zip(fieldnames(typeof(samples)), samples)])

# ╔═╡ 4f97b572-f9b0-11ea-0a99-87af0797bf28
md"""
**Random letters from the alphabet:**
"""

# ╔═╡ 4e8d327e-f9b0-11ea-3f16-c178d96d07d9
md"""
**Random letters at the correct frequencies:**
"""

# ╔═╡ d83f8bbc-f9af-11ea-2392-c90e28e96c65
md"""
**Random letters at the correct transition frequencies:**
"""

# ╔═╡ 0e465160-f937-11ea-0ebb-b7e02d71e8a8
function sample_text(A, n)
	
	first_index = rand_sample(vec(sum(A, dims=1)))
	
	indices = reduce(1:n; init=[first_index]) do word, _
		prev = last(word)
		freq = normalize_array(A[prev, :])
		next = rand_sample(freq)
		[word..., next]
	end
	
	String(alphabet[indices])
end

# ╔═╡ 141af892-f933-11ea-1e5f-154167642809
md"""
It looks like we have a decent language model, in the sense that it understands _transition frequencies_ in the language. In the demo above, try switching the language between $(join(string.(fieldnames(typeof(samples))), " and ")) -- the generated text clearly looks more like one or the other, demonstrating that the model can capture differences between the two languages. What's remarkable is that our "training data" was just a single paragraph per language.

In this exercise, we will use our model to write a **classifier**: a program that automatically classifies a text as either $(join(string.(fieldnames(typeof(samples))), " or ")). 

This is not a difficult task - you can get dictionaries for both languages, and count matches - but we are doing something much more cool: we only use a single paragraph of each language, and we use a _language model_ as classifier.
"""

# ╔═╡ 7eed9dde-f931-11ea-38b0-db6bfcc1b558
html"<h4 id='mystery-detect'>Mystery sample</h4>
<p>Enter some text here -- we will detect in which language it is written!</p>" # dont delete me

# ╔═╡ 7e3282e2-f931-11ea-272f-d90779264456
@bind mystery_sample TextField((70,8), default="""
Small boats are typically found on inland waterways such as rivers and lakes, or in protected coastal areas. However, some boats, such as the whaleboat, were intended for use in an offshore environment. In modern naval terms, a boat is a vessel small enough to be carried aboard a ship. Anomalous definitions exist, as lake freighters 1,000 feet (300 m) long on the Great Lakes are called "boats". 
""")

# ╔═╡ fc58ea34-fedc-11ea-38c1-1dd1d5d6cb53
md""" Original sample

Small boats are typically found on inland waterways such as rivers and lakes, or in protected coastal areas. However, some boats, such as the whaleboat, were intended for use in an offshore environment. In modern naval terms, a boat is a vessel small enough to be carried aboard a ship. Anomalous definitions exist, as lake freighters 1,000 feet (300 m) long on the Great Lakes are called "boats". 
"""

# ╔═╡ 7df55e6c-f931-11ea-33b8-fdc3be0b6cfa
mystery_sample

# ╔═╡ 292e0384-fb57-11ea-0238-0fbe416fc976
md"""
Let's compute the transition frequencies of our mystery sample! Type some text in the box below, and check whether the frequency matrix updates.
"""

# ╔═╡ 7dabee08-f931-11ea-0cb2-c7d5afd21551
transition_frequencies(mystery_sample)

# ╔═╡ 3736a094-fb57-11ea-1d39-e551aae62b1d
md"""
Our model will **compare the transition frequencies of our mystery sample** to those of our two language sample. The closest match will be our detected language.

The only question left is: How do we compare two matrices? When two matrices are almost equal, but not exactly, we want to quantify their _distance_.

👉 Write a function called `matrix_distance` which takes 2 matrices of the same size and finds the distance between them by:

1. Subtracting corresponding elements
2. Finding the absolute value of the difference
3. Summing the differences
"""

# ╔═╡ 13c89272-f934-11ea-07fe-91b5d56dedf8
function matrix_distance(A, B)
	sqrt(sum((A .- B).^2))
	#issing # do something with A .- B
end

# ╔═╡ 7d60f056-f931-11ea-39ae-5fa18a955a77
distances = map(samples) do sample
	matrix_distance(transition_frequencies(mystery_sample), transition_frequencies(sample))
end

# ╔═╡ 7d1439e6-f931-11ea-2dab-41c66a779262
try
	@assert !ismissing(distances.English)
	"""<h2>It looks like this text is **$(argmin(distances))**!</h2>""" |> HTML
catch
end

# ╔═╡ 8c7606f0-fb93-11ea-0c9c-45364892cbb8
md"""
We have written a cell that selects the language with the _smallest distance_ to the mystery language. Make sure sure that `matrix_distance` is working correctly, and [scroll up](#mystery-detect) to the mystery text to see it in action!

#### Further reading
It turns out that the SVD of the transition matrix can mysteriously group the alphabet into vowels and consonants, without any extra information. See [this paper](http://languagelog.ldc.upenn.edu/myl/Moler1983.pdf) if you want to try it yourself! We found that removing the space from `alphabet` (to match the paper) gave better results.
"""

# ╔═╡ 82e0df62-fb54-11ea-3fff-b16c87a7d45b
md"""
## **Exercise 2** - _Language generation_

Our model from Exercise 1 has the property that it can easily be 'reversed' to _generate_ text. While this is useful to demonstrate its structure, the produced text is mostly meaningless: it fails to model words, let alone sentence structure.

To take our model one step further, we are going to _generalize_ what we have done so far. Instead of looking at _letter combinations_, we will model _word combinations_.  And instead of analyzing the frequencies of bigrams (combinations of two letters), we are going to analyze _$n$-grams_.

#### Dataset
This also means that we are going to need a larger dataset to train our model on: the number of English words (and their combinations) is much higher than the number of letters.

We will train our model on the novel [_Emma_ (1815), by Jane Austen](https://en.wikipedia.org/wiki/Emma_(novel)). This work is in the public domain, which means that we can download the whole book as a text file from `archive.org`. We've done the process of downloading and cleaning already, and we have split the text into word and punctuation tokens.
"""

# ╔═╡ b7601048-fb57-11ea-0754-97dc4e0623a1
emma = let
	raw_text = read(download("https://ia800303.us.archive.org/24/items/EmmaJaneAusten_753/emma_pdf_djvu.txt"), String)
	
	first_words = "Emma Woodhouse"
	last_words = "THE END"
	start_index = findfirst(first_words, raw_text)[1]
	stop_index = findlast(last_words, raw_text)[end]
	
	raw_text[start_index:stop_index]
end;

# ╔═╡ cc42de82-fb5a-11ea-3614-25ef961729ab
function splitwords(text)
	# clean up whitespace
	cleantext = replace(text, r"\s+" => " ")
	
	# split on whitespace or other word boundaries
	tokens = split(cleantext, r"(\s|\b)")
end

# ╔═╡ d66fe2b2-fb5a-11ea-280f-cfb12b8296ac
emma_words = splitwords(emma)

# ╔═╡ 4ca8e04a-fb75-11ea-08cc-2fdef5b31944
forest_words = splitwords(samples.English)

# ╔═╡ 6f613cd2-fb5b-11ea-1669-cbd355677649
md"""
#### Exercise 2.1 - _bigrams revisited_

The goal of the upcoming exercises is to **generalize** what we have done in Exercise 1. To keep things simple, we _split up our problem_ into smaller problems. (The solution to any computational problem.)

First, here is a function that takes an array, and returns the array of all **neighbour pairs** from the original. For example,

```julia
bigrams([1, 2, 3, 42])
```
gives

```julia
[[1,2], [2,3], [3,42]]
```

(We used integers as "words" in this example, but our function works with any type of word.)
"""

# ╔═╡ 91e87974-fb78-11ea-3ce4-5f64e506b9d2
function bigrams(words)
	map(1:length(words)-1) do i
		words[i:i+1]
	end
end

# ╔═╡ 9f98e00e-fb78-11ea-0f6c-01206e7221d6
bigrams([1, 2, 3, 42])

# ╔═╡ d7d8cd0c-fb6a-11ea-12bf-2d1448b38162
md"""
👉 Next, it's your turn to write a more general function `ngrams` that takes an array and a number $n$, and returns all **subsequences of length $n$**. For example:

```julia
ngrams([1, 2, 3, 42], 3)
```
should give

```julia
[[1,2,3], [2,3,42]]
```

and

```julia
ngrams([1, 2, 3, 42], 2) == bigrams([1, 2, 3, 42])
```
"""

# ╔═╡ 7be98e04-fb6b-11ea-111d-51c48f39a4e9
function ngrams(words, n)
	map(1:length(words)-n+1) do i
		words[i:i+n-1]
	end
end

# ╔═╡ 052f822c-fb7b-11ea-382f-af4d6c2b4fdb
ngrams([1, 2, 3, 42], 3)

# ╔═╡ 067f33fc-fb7b-11ea-352e-956c8727c79f
ngrams(forest_words, 4)

# ╔═╡ 7b10f074-fb7c-11ea-20f0-034ddff41bc3
md"""

#### Exercise 2.2 - _frequency matrix revisisted_
In Exercise 1, we use a 2D array to store the bigram frequencies, where each column or row corresponds to a character from the alphabet. If we use trigrams, we could store the frequencies in a 3D array, and so on. 

However, when counting words instead of letters, we run into a problem. A 3D array with one row, column and layer per word has too many elements to store on our computer.
"""

# ╔═╡ 24ae5da0-fb7e-11ea-3480-8bb7b649abd5
md"""
_Emma_ consists of $(
	length(Set(emma_words))
) unique words. This means that there are $(
	Int(floor(length(Set(emma_words))^3 / 10^9))
) billion possible trigrams - that's too much!
"""

# ╔═╡ 47836744-fb7e-11ea-2305-3fa5819dc154
md"""
$(html"<br>")

Although the frequency array would be very large, most entries are zero. For example, _"Emma"_ is a common word, but _"Emma Emma Emma"_ does not occur in the novel. This _sparsity_ of non-zero entries can be used to **store the same information in a more efficient structure**. 

Julia's [`SparseArrays.jl` package](https://docs.julialang.org/en/v1/stdlib/SparseArrays/index.html) might sound like a logical choice, but these arrays only support 1D and 2D types, and we also want to directly index using strings, not just integers. So instead, we will use `Dict`, the dictionary type.
"""

# ╔═╡ df4fc31c-fb81-11ea-37b3-db282b36f5ef
healthy = Dict("fruits" => ["🍎", "🍊"], "vegetables" => ["🌽", "🎃", "🍕"])

# ╔═╡ c83b1770-fb82-11ea-20a6-3d3a09606c62
healthy["fruits"]

# ╔═╡ 52970ac4-fb82-11ea-3040-8bd0590348d2
md"""
(Did you notice something funny? The dictionary is _unordered_, this is why the entries were printed in reverse from the definition.)

You can dynamically add or change values of a `Dict` by assigning to `my_dict[key]`. You can check whether a key already exists using `haskey(my_dict, key)`.

👉 Use these two techniques to write a function `word_counts` that takes an array of words, and returns a `Dict` with entries `word => number_of_occurences`.

For example:
```julia
word_counts(["to", "be", "or", "not", "to", "be"])
```
should return
```julia
Dict(
	"to" => 2, 
	"be" => 2, 
	"or" => 1, 
	"not" => 1
)
```
"""

# ╔═╡ 8ce3b312-fb82-11ea-200c-8d5b12f03eea
function word_counts(words::Vector)
	counts = Dict()
	
	for x in words
		if haskey(counts, x)
			counts[x] += 1
		else
			counts[x] = 1
		end
	end
	return counts
end

# ╔═╡ a2214e50-fb83-11ea-3580-210f12d44182
word_counts(["to", "be", "or", "not", "to", "be"])

# ╔═╡ 808abf6e-fb84-11ea-0785-2fc3f1c4a09f
md"""
How many times does `"Emma"` occur in the book?
"""

# ╔═╡ 953363dc-fb84-11ea-1128-ebdfaf5160ee
emma_count = word_counts(emma_words)["Emma"]

# ╔═╡ 294b6f50-fb84-11ea-1382-03e9ab029a2d
md"""
Great! Let's get back to our ngrams. For the purpose of generating text, we are going to store a _completions cache_. This is a dictionary where the keys are $(n-1)$-grams, and the values are all found words that complete it to an $n$-gram. Let's look at an example:

```julia
let
	trigrams = ngrams(split("to be or not to be that is the question", " "), 3)
	cache = completions_cache(trigrams)
	cache == Dict(
		["to", "be"] => ["or", "that"],
		["be", "or"] => ["not"],
		["or", "not"] => ["to"],
		...
	)
end
```

So for trigrams, our keys are the first $2$ words of each trigram, and the values are arrays containing every third word of those trigrams.

If the same ngram occurs multiple times (e.g. "said Emma laughing"), then the last word ("laughing") should also be stored multiple times. This will allow us to generate trigrams with the same frequencies as the original text.

👉 Write the function `completions_cache`, which takes an array of ngrams (i.e. an array of arrays of words, like the result of your `ngram` function), and returns a dictionary like described above.
"""

# ╔═╡ b726f824-fb5e-11ea-328e-03a30544037f
function completions_cache(grams)
	cache = Dict()
	
	for x in grams
		gram_prefix = x[1:lastindex(x)-1]
		gram_suffix = last(x)
		if haskey(cache, gram_prefix)
			push!(cache[gram_prefix],gram_suffix)
		else
			cache[gram_prefix] = [gram_suffix]
		end
	end
	return cache
end

# ╔═╡ 18355314-fb86-11ea-0738-3544e2e3e816
let
	trigrams = ngrams(split("to be or not to be that is the question", " "), 3)
	completions_cache(trigrams)
end

# ╔═╡ 3d105742-fb8d-11ea-09b0-cd2e77efd15c
md"""
#### Exercise 2.4 - _write a novel_

We have everything we need to generate our own novel! The final step is to sample random ngrams, in a way that each next ngram overlaps with the previous one. We've done this in the function `generate_from_ngrams` below - feel free to look through the code, or to implement your own version.
"""

# ╔═╡ a72fcf5a-fb62-11ea-1dcc-11451d23c085
"""
	generate_from_ngrams(grams, num_words)

Given an array of ngrams (i.e. an array of arrays of words), generate a sequence of `num_words` words by sampling random ngrams.
"""
function generate_from_ngrams(grams, num_words)
	n = length(first(grams))
	cache = completions_cache(grams)
	
	# we need to start the sequence with at least n-1 words.
	# a simple way to do so is to pick a random ngram!
	sequence = [rand(grams)...]
	
	# we iteratively add one more word at a time
	for i ∈ n+1:num_words
		# the previous n-1 words
		tail = sequence[end-(n-2):end]
		
		# possible next words
		completions = cache[tail]
		
		choice = rand(completions)
		push!(sequence, choice)
	end
	sequence
end

# ╔═╡ f83991c0-fb7c-11ea-0e6f-1f80709d00c1
"Compute the ngrams of an array of words, but add the first n-1 at the end, to ensure that every ngram ends in the the beginning of another ngram."
function ngrams_circular(words, n)
	ngrams([words..., words[1:n-1]...], n)
end

# ╔═╡ abe2b862-fb69-11ea-08d9-ebd4ba3437d5
completions_cache(ngrams_circular(forest_words, 3))

# ╔═╡ 4b27a89a-fb8d-11ea-010b-671eba69364e
"""
	generate(source_text::AbstractString, num_token; n=3, use_words=true)

Given a source text, generate a `String` that "looks like" the original text by satisfying the same ngram frequency distribution as the original.
"""
function generate(source_text::AbstractString, s; n=3, use_words=true)
	preprocess = if use_words
		splitwords
	else
		collect
	end
	
	words = preprocess(source_text)
	if length(words) < n
		""
	else
		grams = ngrams_circular(words, n)
		result = generate_from_ngrams(grams, s)
		if use_words
			join(result, " ")
		else
			String(result)
		end
	end
end

# ╔═╡ d7b7a14a-fb90-11ea-3e2b-2fd8f379b4d8
md"
#### Interactive demo

Enter your own text in the box below, and use that as training data to generate anything!
"

# ╔═╡ 1939dbea-fb63-11ea-0bc2-2d06b2d4b26c
@bind generate_demo_sample TextField((50,5), default=samples.English)

# ╔═╡ 70169682-fb8c-11ea-27c0-2dad2ff3080f
md"""Using $(@bind generate_sample_n_letters NumberField(1:5))grams for characters"""

# ╔═╡ 402562b0-fb63-11ea-0769-375572cc47a8
md"""Using $(@bind generate_sample_n_words NumberField(1:5))grams for words"""

# ╔═╡ 2521bac8-fb8f-11ea-04a4-0b077d77529e
md"""
### Automatic Jane Austen

Uncomment the cell below to generate some Jane Austen text:
"""

# ╔═╡ 6b4d6584-f3be-11ea-131d-e5bdefcc791b
md"## Function library

Just some helper functions used in the notebook."

# ╔═╡ 54b1e236-fb53-11ea-3769-b382ef8b25d6
function Quote(text::AbstractString)
	text |> Markdown.Paragraph |> Markdown.BlockQuote |> Markdown.MD
end

# ╔═╡ be55507c-f9a7-11ea-189c-4ffe8377212e
if sample_freqs !== missing
	String([rand_sample_letter(sample_freqs) for _ in 1:400]) |> Quote
end

# ╔═╡ 46c905d8-f9b0-11ea-36ed-0515e8ed2621
String(rand(alphabet, 400)) |> Quote

# ╔═╡ 489b03d4-f9b0-11ea-1de0-11d4fe4e7c69
String([rand_sample_letter(letter_frequencies(ex23_sample)) for _ in 1:400]) |> Quote

# ╔═╡ fd202410-f936-11ea-1ad6-b3629556b3e0
sample_text(transition_frequencies(clean(ex23_sample)), 400) |> Quote

# ╔═╡ b5dff8b8-fb6c-11ea-10fc-37d2a9adae8c
generate(
	generate_demo_sample, 400; 
	n=generate_sample_n_letters, 
	use_words=false
) |> Quote

# ╔═╡ ee8c5808-fb5f-11ea-19a1-3d58217f34dc
generate(
	generate_demo_sample, 100; 
	n=generate_sample_n_words, 
	use_words=true
) |> Quote

# ╔═╡ 49b69dc2-fb8f-11ea-39af-030b5c5053c3
 generate(emma, 100; n=4) |> Quote

# ╔═╡ ddef9c94-fb96-11ea-1f17-f173a4ff4d89
function compimg(img, labels=[c*d for c in replace(alphabet, ' ' => "_"), d in replace(alphabet, ' ' => "_")])
	xmax, ymax = size(img)
	xmin, ymin = 0, 0
	arr = [(j-1, i-1) for i=1:ymax, j=1:xmax]

	compose(context(units=UnitBox(xmin, ymin, xmax, ymax)),
		fill(vec(img)),
		compose(context(),
			fill("white"), font("monospace"), 
			text(first.(arr) .+ .1, last.(arr) .+ 0.6, labels)),
		rectangle(
			first.(arr),
			last.(arr),
			fill(1.0, length(arr)),
			fill(1.0, length(arr))))
end

# ╔═╡ b7803a28-fb96-11ea-3e30-d98eb322d19a
function show_pair_frequencies(A)
	imshow = let
		to_rgb(x) = RGB(0.36x, 0.82x, 0.8x)
		to_rgb.(A ./ maximum(abs.(A)))
	end
	compimg(imshow)
end

# ╔═╡ ace3dc76-f9ae-11ea-2bee-3d0bfa57cfbc
show_pair_frequencies(transition_frequencies(first_sample))

# ╔═╡ 00115b6e-f381-11ea-0bc6-61ca119cb628
bigbreak = html"<br><br><br><br><br>";

# ╔═╡ c086bd1e-f384-11ea-3b26-2da9e24360ca
bigbreak

# ╔═╡ eaa8c79e-f9a2-11ea-323f-8bb2bd36e11c
md"""
$(bigbreak)
#### Exercise 1.2 - _Letter frequencies_

We are going to count the _frequency_ of each letter in this sample, after applying your `clean` function. Can you guess which character is most frequent?
"""

# ╔═╡ dcffd7d2-f9a6-11ea-2230-b1afaecfdd54
md"""
$(bigbreak)
Now that we know the frequencies of letters in English, we can generate random text that already looks closer to English!

**Random letters from the alphabet:**
"""

# ╔═╡ 77623f3e-f9a9-11ea-2f46-ff07bd27cd5f
md"""
$(bigbreak)
#### Exercise 1.3 - _Transition frequencies_
In the previous exercise we computed the frequency of each letter in the sample by _counting_ their occurrences, and then dividing by the total number of counts.

In this exercise, we are going to count _letter transitions_, such as `aa`, `as`, `rt`, `yy`. Two letters might both be common, like `a` and `e`, but their combination, `ae`, is uncommon in English. 

To quantify this observation, we will do the same as in our last exercise: we count occurrences in a _sample text_, to create the **transition frequency matrix**.
"""

# ╔═╡ d3d7bd9c-f9af-11ea-1570-75856615eb5d
bigbreak

# ╔═╡ 6718d26c-f9b0-11ea-1f5a-0f22f7ddffe9
md"""
$(bigbreak)

#### Exercise 1.4 - _Language detection_
"""

# ╔═╡ 568f0d3a-fb54-11ea-0f77-171718ef12a5
bigbreak

# ╔═╡ 7f341c4e-fb54-11ea-1919-d5421d7a2c75
bigbreak

# ╔═╡ Cell order:
# ╠═85cfbd10-f384-11ea-31dc-b5693630a4c5
# ╟─938185ec-f384-11ea-21dc-b56b7469f798
# ╠═86e1ee96-f314-11ea-03f6-0f549b79e7c9
# ╠═a4937996-f314-11ea-2ff9-615c888afaa8
# ╟─c086bd1e-f384-11ea-3b26-2da9e24360ca
# ╟─3cc688d2-f996-11ea-2a6f-0b4c7a5b74c2
# ╠═829eea8e-05f0-11eb-2a9f-cf1b0c339053
# ╠═28195e06-05f4-11eb-3997-2377c3e3b124
# ╠═9424c438-0652-11eb-2329-4b90f2f7cecf
# ╠═ae61d73c-0652-11eb-38b5-df004da5eb69
# ╠═b6e8e120-0652-11eb-3af6-dbc37e3ae9b5
# ╠═0d3aaa6a-0651-11eb-0034-d77db75bb34c
# ╠═d67034d0-f92d-11ea-31c2-f7a38ebb412f
# ╟─a094e2ac-f92d-11ea-141a-3566552dd839
# ╟─faf80910-0650-11eb-0798-dd4b1628700e
# ╟─5c74a052-f92e-11ea-2c5b-0f1a3a14e313
# ╟─05f0182c-f999-11ea-0a52-3d46c65a049e
# ╟─24860970-fc48-11ea-0009-cddee695772c
# ╟─734851c6-f92d-11ea-130d-bf2a69e89255
# ╟─8d3bc9ea-f9a1-11ea-1508-8da4b7674629
# ╠═4affa858-f92e-11ea-3ece-258897c37e51
# ╟─eaa8c79e-f9a2-11ea-323f-8bb2bd36e11c
# ╠═2680b506-f9a3-11ea-0849-3989de27dd9f
# ╟─571d28d6-f960-11ea-1b2e-d5977ecbbb11
# ╠═6a64ab12-f960-11ea-0d92-5b88943cdb1a
# ╟─603741c2-f9a4-11ea-37ce-1b36ecc83f45
# ╟─b3de6260-f9a4-11ea-1bae-9153a92c3fe5
# ╠═a6c36bd6-f9a4-11ea-1aba-f75cecc90320
# ╟─6d3f9dae-f9a5-11ea-3228-d147435e266d
# ╠═92bf9fd2-f9a5-11ea-25c7-5966e44db6c6
# ╠═dcffd7d2-f9a6-11ea-2230-b1afaecfdd54
# ╟─01215e9a-f9a9-11ea-363b-67392741c8d4
# ╟─be55507c-f9a7-11ea-189c-4ffe8377212e
# ╟─8ae13cf0-f9a8-11ea-3919-a735c4ed9e7f
# ╟─343d63c2-fb58-11ea-0cce-efe1afe070c2
# ╟─b5b8dd18-f938-11ea-157b-53b145357fd1
# ╟─0e872a6c-f937-11ea-125e-37958713a495
# ╟─77623f3e-f9a9-11ea-2f46-ff07bd27cd5f
# ╠═fbb7c04e-f92d-11ea-0b81-0be20da242c8
# ╠═80118bf8-f931-11ea-34f3-b7828113ffd8
# ╠═7f4f6ce4-f931-11ea-15a4-b3bec6a7e8b6
# ╠═d40034f6-f9ab-11ea-3f65-7ffd1256ae9d
# ╟─689ed82a-f9ae-11ea-159c-331ff6660a75
# ╠═ace3dc76-f9ae-11ea-2bee-3d0bfa57cfbc
# ╟─0b67789c-f931-11ea-113c-35e5edafcbbf
# ╠═6896fef8-f9af-11ea-0065-816a70ba9670
# ╟─39152104-fc49-11ea-04dd-bb34e3600f2f
# ╟─d3d7bd9c-f9af-11ea-1570-75856615eb5d
# ╟─2f8dedfc-fb98-11ea-23d7-2159bdb6a299
# ╟─b7446f34-f9b1-11ea-0f39-a3c17ba740e5
# ╟─4f97b572-f9b0-11ea-0a99-87af0797bf28
# ╟─46c905d8-f9b0-11ea-36ed-0515e8ed2621
# ╟─4e8d327e-f9b0-11ea-3f16-c178d96d07d9
# ╟─489b03d4-f9b0-11ea-1de0-11d4fe4e7c69
# ╟─d83f8bbc-f9af-11ea-2392-c90e28e96c65
# ╟─fd202410-f936-11ea-1ad6-b3629556b3e0
# ╟─0e465160-f937-11ea-0ebb-b7e02d71e8a8
# ╟─6718d26c-f9b0-11ea-1f5a-0f22f7ddffe9
# ╟─141af892-f933-11ea-1e5f-154167642809
# ╟─7eed9dde-f931-11ea-38b0-db6bfcc1b558
# ╟─7e3282e2-f931-11ea-272f-d90779264456
# ╠═fc58ea34-fedc-11ea-38c1-1dd1d5d6cb53
# ╟─7d1439e6-f931-11ea-2dab-41c66a779262
# ╠═7df55e6c-f931-11ea-33b8-fdc3be0b6cfa
# ╟─292e0384-fb57-11ea-0238-0fbe416fc976
# ╠═7dabee08-f931-11ea-0cb2-c7d5afd21551
# ╟─3736a094-fb57-11ea-1d39-e551aae62b1d
# ╠═13c89272-f934-11ea-07fe-91b5d56dedf8
# ╟─7d60f056-f931-11ea-39ae-5fa18a955a77
# ╟─8c7606f0-fb93-11ea-0c9c-45364892cbb8
# ╟─568f0d3a-fb54-11ea-0f77-171718ef12a5
# ╟─82e0df62-fb54-11ea-3fff-b16c87a7d45b
# ╠═b7601048-fb57-11ea-0754-97dc4e0623a1
# ╟─cc42de82-fb5a-11ea-3614-25ef961729ab
# ╠═d66fe2b2-fb5a-11ea-280f-cfb12b8296ac
# ╠═4ca8e04a-fb75-11ea-08cc-2fdef5b31944
# ╟─6f613cd2-fb5b-11ea-1669-cbd355677649
# ╠═91e87974-fb78-11ea-3ce4-5f64e506b9d2
# ╠═9f98e00e-fb78-11ea-0f6c-01206e7221d6
# ╟─d7d8cd0c-fb6a-11ea-12bf-2d1448b38162
# ╠═7be98e04-fb6b-11ea-111d-51c48f39a4e9
# ╠═052f822c-fb7b-11ea-382f-af4d6c2b4fdb
# ╠═067f33fc-fb7b-11ea-352e-956c8727c79f
# ╠═7b10f074-fb7c-11ea-20f0-034ddff41bc3
# ╟─24ae5da0-fb7e-11ea-3480-8bb7b649abd5
# ╟─47836744-fb7e-11ea-2305-3fa5819dc154
# ╠═df4fc31c-fb81-11ea-37b3-db282b36f5ef
# ╠═c83b1770-fb82-11ea-20a6-3d3a09606c62
# ╟─52970ac4-fb82-11ea-3040-8bd0590348d2
# ╠═8ce3b312-fb82-11ea-200c-8d5b12f03eea
# ╠═a2214e50-fb83-11ea-3580-210f12d44182
# ╟─808abf6e-fb84-11ea-0785-2fc3f1c4a09f
# ╠═953363dc-fb84-11ea-1128-ebdfaf5160ee
# ╟─294b6f50-fb84-11ea-1382-03e9ab029a2d
# ╠═b726f824-fb5e-11ea-328e-03a30544037f
# ╠═18355314-fb86-11ea-0738-3544e2e3e816
# ╠═abe2b862-fb69-11ea-08d9-ebd4ba3437d5
# ╟─3d105742-fb8d-11ea-09b0-cd2e77efd15c
# ╟─a72fcf5a-fb62-11ea-1dcc-11451d23c085
# ╟─f83991c0-fb7c-11ea-0e6f-1f80709d00c1
# ╟─4b27a89a-fb8d-11ea-010b-671eba69364e
# ╟─d7b7a14a-fb90-11ea-3e2b-2fd8f379b4d8
# ╟─1939dbea-fb63-11ea-0bc2-2d06b2d4b26c
# ╟─70169682-fb8c-11ea-27c0-2dad2ff3080f
# ╟─b5dff8b8-fb6c-11ea-10fc-37d2a9adae8c
# ╟─402562b0-fb63-11ea-0769-375572cc47a8
# ╟─ee8c5808-fb5f-11ea-19a1-3d58217f34dc
# ╟─2521bac8-fb8f-11ea-04a4-0b077d77529e
# ╠═49b69dc2-fb8f-11ea-39af-030b5c5053c3
# ╟─7f341c4e-fb54-11ea-1919-d5421d7a2c75
# ╟─6b4d6584-f3be-11ea-131d-e5bdefcc791b
# ╟─54b1e236-fb53-11ea-3769-b382ef8b25d6
# ╟─b7803a28-fb96-11ea-3e30-d98eb322d19a
# ╟─ddef9c94-fb96-11ea-1f17-f173a4ff4d89
# ╟─00115b6e-f381-11ea-0bc6-61ca119cb628
