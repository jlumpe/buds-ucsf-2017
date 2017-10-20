import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.plotly as py


def get_fitness_by_effect(fitness, effect):
    fitness["position"] = list(range(140))
    fitness_score_list = []
    fitness_score_res = []

    for i in range(len(effect)):
        wt, codon, mut = get_mutation(str(effect.iloc[i]["Variant"]))
        fitness_score_list.append(fitness.iloc[int(codon)-1][mut])
        fitness_score_res.append(codon)
    return fitness_score_list, fitness_score_res


def get_mutation(name):
    return name[0], name[1:-1], name[-1]


def get_effect(fitness, quality, category):
    return fitness[fitness[quality] == category]


def get_average_fit(fitness):
    fitness["average"] = fitness.mean(axis=1)
    return fitness

if __name__ == "__main__":

    fitness = pd.read_csv("fitness.csv")
    effect_prediction = pd.read_csv("wt_asyn_PP_csv.csv")

    get_average_fit(fitness)

    effect = get_effect(effect_prediction, "Predicted Effect", "effect")
    no_effect = get_effect(effect_prediction, "Predicted Effect", "neutral")

    effect_fitness, effect_res = get_fitness_by_effect(fitness, effect)
    no_effect_fitness, no_effect_res = get_fitness_by_effect(fitness, no_effect)

    #######################################################################

    # scatterplot of mean fitness by residue
    # using matplotlib
    # fig, ax = plt.subplots()
    # plt.style.use('fivethirtyeight')
    # points = ax.scatter(fitness.index, fitness["average"], c=fitness.index, s=20, cmap="bwr")
    # fig.colorbar(points, ax=ax)
    # plt.xlabel("Residue Number")
    # plt.ylabel("Mean Fitness")
    # # plt.show()
    # plt.savefig("scatterplot.png")

    #######################################################################
    # bar plot of mean fitness by effect of mutation
    # using plotly

    # code used to communicate with API
    # trace0 = go.Box(
    #     x=effect_fitness,
    #     name='Effect',
    #     marker=dict(color='rgb(10, 140, 208)',
    # ),
    #     boxmean=True
    # )
    #
    # trace1 = go.Box(
    #     x=no_effect_fitness,
    #     name='No Effect',
    #     marker=dict(color='rgb(10, 140, 208)',
    #                 ),
    #     boxmean=True
    # )

    # data = [trace0, trace1]
    # py.iplot(data)

    # # to use the final version and save
    # fig = py.get_figure("https://plot.ly/~christacaggiano/10/")
    # py.image.save_as(fig, "boxplot.png")

    #######################################################################

    fig, ax = plt.subplots()

    effect_fitness = pd.DataFrame({"fit":effect_fitness, "res":effect_res, "effect":"effect"})
    no_effect_fitness=pd.DataFrame({"fit":no_effect_fitness, "res":no_effect_res, "effect":"no_effect"})
    concat = pd.concat([effect_fitness, no_effect_fitness])

    sns.boxplot(y="fit", x="effect", data=concat, palette=sns.light_palette("grey"))
    sns.stripplot(y="fit", x="effect", data=concat, hue="res", jitter=True, palette=sns.color_palette("RdBu_r", 140))
    ax.legend_.remove()
    plt.show()


