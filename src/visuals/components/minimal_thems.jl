function theme_minimal_2()
    Theme(
        Axis = (
            backgroundcolor = :transparent,
            xgridvisible = false,
            ygridvisible = false,
            xminorgridvisible = false,
            yminorgridvisible = false,
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xminorticksvisible = true,
            yminorticksvisible = true,
            xticksvisible = false,
            yticksvisible = false,
            xlabelpadding = 3,
            ylabelpadding = 3
        ),
        Legend = (
            framevisible = false,
            padding = (0, 0, 0, 0),
        ),
        Axis3 = (
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xspinesvisible = true,
            yspinesvisible = true,
            zspinesvisible = true,
            yzpanelvisible = false,
            xzpanelvisible = false,
            xypanelvisible = false,
            xticksvisible = true,
            yticksvisible = true,
            zticksvisible = false,
            xticklabelpad = 3,
            yticklabelpad = 3,
            zticklabelpad = 6,
            xspinecolor_2 = :transparent,
            xspinecolor_3 = :transparent,
            yspinecolor_2 = :transparent,
            yspinecolor_3 = :transparent,
            zspinecolor_2 = :transparent,
            zspinecolor_3 = :transparent,
        ),
        Colorbar = (
            ticksvisible = false,
            spinewidth = 0,
            ticklabelpad = 5,
        )
    )
end