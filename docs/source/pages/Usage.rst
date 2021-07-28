Usage
=========================================================

A simple plot using `imshow` can be obtained by::


        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.imshow(mean_curvature, cmap='bwr', interpolation='gaussian', origin='lower')
        ax.set_title('Mean Curvature')
        plt.show()


As an alternative, you can use contour plots using `contourf`::


        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.contourf(mean_curvature, cmap='bwr, origin='lower')
        ax.set_title('Mean Curvature')
        plt.show()
