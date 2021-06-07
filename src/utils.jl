import Printf: @sprintf

"""Returns a human readable string representing seconds in terms of days, hours, minutes or seconds."""
function readable_secs(secs::Real)
    days = Int(floor(secs/3600/24))
    hours = Int(floor((secs/3600) % 24))
    minutes = Int(floor((secs/60) % 60))
    seconds = Int(floor(secs%3600%60))
    secs1f = @sprintf "%.1fs" secs%3600%60
    secs2f = @sprintf "%.2fs" secs%3600%60

    if days > 0
        return "$(days)d, $(hours)h"
    elseif hours > 0
        return "$(hours)h, $(minutes)min"
    elseif minutes > 0
        return "$(minutes)min, $(seconds)s"
    elseif seconds > 10
        return secs1f
    else
        return secs2f
    end
end