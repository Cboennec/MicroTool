
# Load packages
include("..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames


# My model
x = 1.0:-0.01:0.0
recycle_value = 10
full_price = 150
EOL = 0.8


# Décroissance exponentielle pour arriver à environ 20% de valeur utile restante à EOL. on y ajout le prix de recyclage
# Les deux formules suivantes sont des propositions appuyer sur des arguments logique mais calibré sur aucune donnée ni savoir.

# Le principe est le suivant :
# La batterie pleine vaux 100% de sa valeur, à EOL on la change, on considère alors qu'elle va bientôt être inutile,
# Néanmoins même inutilisable elle conserve son prix de recyclage des matériaux.
# Mon pref  c'est le 1, le plus raisonnable c'est probablement le 3.

value =  x .^ (2 ./ (1 - EOL))  .* full_price .*  (1 - (recycle_value/full_price))   .+    recycle_value
value2 = exp.(2 * (x .- 1) ./ (1 - EOL))  .* full_price .*  (1 - (recycle_value/full_price))   .+    recycle_value
value3 = [(1 - ((1-a)/(1-EOL))) >= 0 ? (1 - ((1-a)/(1-EOL))) : 0 for a in x]  .* full_price .*  (1 - (recycle_value/full_price)) .+ recycle_value

model_names = [repeat(["x^y"], length(value)); repeat(["exp(x)"], length(value2)); repeat(["linear"],length(value3)) ]
X = repeat(x,3)

values = [value; value2; value3]
data = Pandas.DataFrame(Dict("value" => values, "x" => X, "Models" => model_names))
myplt = Seaborn.lineplot(data = data, y = "value", x = "x", hue = "Models")

Seaborn.plt.ylim(0,full_price)
Seaborn.plt.xlim(1,0)

myplt.set_xlabel("SOH")
myplt.set_ylabel("Price (€)")
myplt.axvline(x=EOL, ymin=0, ymax=full_price, color="red")
myplt.set_xticks([1.0:-0.1:0.0; EOL])
myplt.set_xticklabels([[string(x) for x in 1.0:-0.1:0.0]; "EOL"])
myplt.axhline(y=recycle_value, xmin=0, xmax=1, color="purple")
myplt.set_yticks([0, recycle_value, full_price])
myplt.set_yticklabels(["0", "Recycle value", "Full price"])
Seaborn.legend(title ="Models")
