
plt.figure(figsize=(10, 4))
ax = plt.gca()
hierarchy.dendrogram(z, no_labels=True, color_threshold=0., above_threshold_color='k')
ax.axis('off')
plt.tight_layout()
plt.savefig('methyl_dendrogram_%dcuts.png' % 0, dpi=200)

for i in range(2, 7):
    c = clustering.dendrogram_threshold_by_nclust(z, i)
    plt.figure(figsize=(10, 4))
    ax = plt.gca()
    hierarchy.dendrogram(z, no_labels=True, color_threshold=c, above_threshold_color='k')
    ax.axhline(c, c='gray', ls='--')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('methyl_dendrogram_%dcuts.png' % i, dpi=200)
