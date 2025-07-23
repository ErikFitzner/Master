using FlameGraphs, FileIO, ProfileView

data, lidict = load("C:/Users/User/Downloads/Profile_DSF.jlprof")   # tuple zurück
g = flamegraph(data; lidict)                   # Flame-Graph bauen
ProfileView.view(g)                            # interaktiv inspizieren
