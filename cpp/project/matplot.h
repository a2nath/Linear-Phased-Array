#pragma once
#include <vector>
#include "common.h"
#include "matplotlibcpp.h"

/* plot the charts */
namespace plt = matplotlibcpp;

void plot(const std::vector<Placements>& input, std::string title = "", bool hold = false)
{

    //#ifndef _DEBUG
        /* plot only */
    std::vector<int> locx, locy;
    for (auto& loc : input)
    {
        locx.emplace_back(loc.x);
        locy.emplace_back(loc.y);
    }

    //if (hold) plt::hold(true);

    plt::scatter(locx, locy, { {"color", "green"}, {"size", "20"} });
    //plt::
    if (title.size()) plt::title(title);
    plt::legend();
    plt::show();
    //#endif
}
