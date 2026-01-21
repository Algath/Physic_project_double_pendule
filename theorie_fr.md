# Pendule double — Théorie et modélisation

*Ce document vise à fournir un support théorique au niveau Bachelor, en lien direct avec l’implémentation numérique du projet.*

Ce document fournit une présentation centrée sur la théorie du pendule double planaire utilisé dans le projet, incluant la cinématique, les expressions d'énergie, une dérivation lagrangienne des équations du mouvement, les approximations courantes (petits angles), les options d'amortissement, les considérations numériques, la stratégie d'estimation des paramètres, et une brève comparaison avec le mouvement balistique.

**Objectif :** fournir des équations claires et reproductibles qui correspondent à l'implémentation dans `pendule.jl` et l'arrière‑plan nécessaire pour calibrer et étendre le modèle numérique.

---

## Qu'est‑ce qu'un pendule double ?

Un pendule double est un système mécanique simple composé de deux pendules attachés l'un à l'autre : un pendule (masse ponctuelle $m_1$ sur une tige de longueur $L_1$) avec un second pendule (masse $m_2$ sur une tige de longueur $L_2$) fixé à son extrémité. On le nomme aussi « pendule composé » ou — quand on insiste sur son comportement dynamique — « pendule chaotique » ou « pendule double planaire ».

Propriétés clés :

- Déterministe mais souvent chaotique : le mouvement est entièrement déterminé par les lois de Newton (ou équivalemment par une formulation lagrangienne), mais pour de nombreuses conditions initiales le système montre une sensibilité aux conditions initiales — de petites variations d'angles ou de vitesses entraînent des trajectoires qui divergent rapidement.
- Faible dimension mais comportement riche : le pendule double planaire minimal a quatre variables d'état (deux angles et deux vitesses angulaires) et peut produire des trajectoires périodiques, quasi‑périodiques ou chaotiques selon l'énergie et les conditions initiales.
- Usages courants : démonstration du chaos et de la dynamique non linéaire, banc d'essai pour intégrateurs numériques, exemples pédagogiques en mécanique. Il apparaît aussi en robotique et contrôle lorsqu'on analyse des segments articulés.

Conséquence pratique pour ce projet : à cause de cette sensibilité, l'estimation des paramètres et les prédictions à court terme doivent être traitées avec précaution (fenêtres d'ajustement courtes, filtrage du bruit, quantification des incertitudes via des ensembles).

## 1. Définition du système et cinématique

Géométrie et notation :

- Deux tiges rigides, de masse négligeable, de longueurs $L_1,L_2$.
- Deux masses ponctuelles $m_1$ (supérieure) et $m_2$ (inférieure).
- Angles $	\theta_1(t),	\theta_2(t)$ mesurés depuis la verticale descendante (donc $	\theta=0$ est la configuration verticale de repos).

Positions cartésiennes (utiles pour relier les données de suivi aux angles) :

$$
\mathbf{r}_1=(x_1,y_1)=(L_1\sin\theta_1,\;L_1\cos\theta_1),
\qquad
\mathbf{r}_2=(x_2,y_2)=(L_1\sin\theta_1+L_2\sin\theta_2,\;L_1\cos\theta_1+L_2\cos\theta_2).
$$

Vitesses (dérivées temporelles) :

$$
\dot{\mathbf{r}}_1=(L_1\cos\theta_1\,\dot\theta_1,\;-L_1\sin\theta_1\,\dot\theta_1),\qquad
\dot{\mathbf{r}}_2=(L_1\cos\theta_1\,\dot\theta_1+L_2\cos\theta_2\,\dot\theta_2,\;-L_1\sin\theta_1\,\dot\theta_1-L_2\sin\theta_2\,\dot\theta_2).
$$

Ces expressions servent à calculer l'énergie cinétique et à transformer les coordonnées de suivi (pixels) en variables angulaires utilisées par le modèle.

## 2. Énergies et Lagrangien

*Les forces de tension dans les tiges sont des forces de contrainte et ne sont pas écrites explicitement dans le formalisme lagrangien.*

- Énergie cinétique :

$$
T=\tfrac12 m_1\lVert\dot{\mathbf{r}}_1\rVert^2+\tfrac12 m_2\lVert\dot{\mathbf{r}}_2\rVert^2.
$$

En développant les normes au carré, on obtient des termes explicites en $\dot\theta_1^2$, $\dot\theta_2^2$ et le couplage $\dot\theta_1\dot\theta_2\cos(\theta_1-\theta_2)$.

- Énergie potentielle (gravitationnelle, référence à l'origine) :

$$
V=m_1 g y_1 + m_2 g y_2 = (m_1+m_2)gL_1\cos\theta_1 + m_2 g L_2\cos\theta_2.
$$

- Lagrangien : $\mathcal{L} = T - V$.

On dérive les équations d'Euler–Lagrange pour les coordonnées généralisées $\theta_1,\theta_2$ :

$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_i}\right)-\frac{\partial\mathcal{L}}{\partial\theta_i}=Q_i,\qquad i=1,2
$$

où $Q_i$ sont les couples généralisés non conservatifs (nuls pour le modèle conservatif).

### Principe d'Hamilton et dérivation Euler–Lagrange

Les équations ci‑dessus découlent du principe d'Hamilton (principe de l'action stationnaire) : la trajectoire réelle $q(t)$ d'un système mécanique entre des temps fixés $t_1$ et $t_2$ rend l'action

$$
S[q]=\int_{t_1}^{t_2}\mathcal{L}(q,\dot q,t)\,dt
$$

stationnaire pour toutes les variations $\delta q(t)$ s'annulant aux extrémités. La variation de l'action est

$$
\delta S=\int_{t_1}^{t_2}\left(\frac{\partial\mathcal{L}}{\partial q}\delta q + \frac{\partial\mathcal{L}}{\partial\dot q}\delta\dot q\right)dt.
$$

En intégrant par parties le second terme et en utilisant $\delta q(t_1)=\delta q(t_2)=0$ on obtient

$$
\delta S=\int_{t_1}^{t_2}\left(\frac{\partial\mathcal{L}}{\partial q}-\frac{d}{dt}\frac{\partial\mathcal{L}}{\partial\dot q}\right)\delta q\,dt.
$$

Comme $\delta q$ est arbitraire, l'intégrande doit s'annuler, d'où l'équation d'Euler–Lagrange pour chaque coordonnée généralisée :

$$
\frac{\partial\mathcal{L}}{\partial q}-\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot q}\right)=0.
$$

Si des forces généralisées non conservatives $Q_i$ agissent, elles apparaissent au membre de droite et donnent la forme utilisée plus haut.

Appliqué à $q=\theta_1,\theta_2$, cette procédure (calculer $T$ et $V$, former $\mathcal{L}=T-V$, calculer les dérivées partielles et temporelles) conduit aux formules explicites pour les accélérations angulaires données en §3. L'approche lagrangienne est pratique car elle prend automatiquement en compte les contraintes cinématiques et utilise directement les coordonnées généralisées.

Le rôle des équations d’Euler–Lagrange n’est pas de résoudre le système analytiquement, mais d’obtenir des équations exploitables numériquement.

## 3. Équations du mouvement explicites (forme finale)

Soit $\delta=\theta_2-\theta_1$, et définissons le dénominateur pratique

$$
\mathrm{den}=m_1+m_2\sin^2\delta.
$$

Alors les accélérations angulaires peuvent s'écrire (à un changement de notation près à l'implémentation dans `pendule.jl`) :

$$
\begin{aligned}
\ddot\theta_1 &= \frac{m_2\sin\delta\big(L_1\dot\theta_1^2\cos\delta + L_2\dot\theta_2^2\big) + m_2 g \sin\theta_2\cos\delta - (m_1 + m_2) g \sin\theta_1}{L_1\,\mathrm{den}},\\[6pt]
\ddot\theta_2 &= \frac{-m_2 L_2\dot\theta_2^2\sin\delta\cos\delta - (m_1 + m_2) L_1\dot\theta_1^2\sin\delta + (m_1 + m_2) g \sin\theta_1\cos\delta - (m_1 + m_2) g \sin\theta_2}{L_2\,\mathrm{den}}.
\end{aligned}
$$

Ces expressions sont non linéaires et couplent angles et vitesses angulaires ; ce sont celles implémentées dans `equations_double_pendulum!` de `pendule.jl`.

### Remarques sur l'algèbre

- Le couplage apparaît via les facteurs $\sin\delta$ et $\cos\delta$ et via les produits des vitesses angulaires au carré.
- Le dénominateur provient de l'élimination des forces de réaction internes et dépend de la géométrie et des masses.

## 4. Approximations et limites courantes

- Linéarisation petits angles (pour des $\theta_i$ petits) : remplacer $\sin\theta\approx\theta$, $\cos\theta\approx1$, ne conserver que les termes linéaires. On obtient un système linéaire 2×2 à coefficients constants ; utile pour analyser les modes normaux et les fréquences propres mais invalide pour de grands mouvements ou les régimes chaotiques.
- Limite pendule simple : poser $m_2\to0$ (ou $L_2\to0$) et récupérer l'équation du pendule simple
  $\ddot\theta + \frac{g}{L}\sin\theta=0.$
- Conservation de l'énergie : dans le modèle conservatif $E=T+V$ reste constant. Une décroissance observée de l'énergie dans les données indique des forces non conservatives (amortissement, frottements), ce qui motive l'ajout de termes d'amortissement et l'estimation de paramètres.

## 5. Amortissement et forces non conservatrices

*Les forces de tension dans les tiges sont des forces de contrainte et ne sont pas écrites explicitement dans le formalisme lagrangien.*

Deux modèles d'amortissement courants :

- Visqueux (couple proportionnel à la vitesse angulaire) : ajouter des couples généralisés $Q_i=-c_i\dot\theta_i$. En pratique, on peut soustraire un terme des formes de $\ddot\theta_i$ de la forme $c_i\dot\theta_i/L_i$ (comme approximation cohérente avec la manière dont les couples se traduisent en accélérations angulaires dans les équations choisies).
- Traînée aérodynamique (plus complexe) : un couple de traînée dépendant de la norme de $\dot\theta$ ou de coefficients fluidodynamiques ; souvent non linéaire en vitesse et éventuellement différent pour chaque masse.

L'ajustement de $c_i$ (ou d'un modèle de perte équivalent) peut se faire en faisant correspondre la décroissance d'énergie observée sur une fenêtre de calibration, ou en optimisant les positions/angles simulés par rapport aux données de suivi sur une fenêtre ultérieure. Commencer par un amortissement visqueux simple est recommandé.

## 6. Intégration numérique et notes d'implémentation

- Le système est un ensemble de 4 EDOs du premier ordre pour $[\theta_1,\omega_1,\theta_2,\omega_2]$ avec le membre de droite donné par les formules ci‑dessus. Utiliser un intégrateur robuste non raide (p.ex. `Tsit5`) depuis `DifferentialEquations.jl` comme dans `pendule.jl`.
- Sauvegarder la solution aux mêmes instants d'échantillonnage que le traceur (`saveat`) pour comparer directement les angles simulés avec les angles observés.
- Divergence chaotique : même avec un modèle parfait, les trajectoires peuvent diverger des mesures à cause de la sensibilité aux conditions initiales. Employer des fenêtres courtes pour l'ajustement des paramètres et évaluer la compétence prédictive de façon probabiliste (ensembles) plutôt que déterministe pour des horizons longs.

## 7. Formulation en variables d'état

Cette relation est avant tout une définition simple et utile :

$$
\dot{\theta}=\omega
$$

Autrement dit, la dérivée de l'angle est la vitesse angulaire. Les autres équations (ci‑dessous notées via la fonction $f$) donnent les accélérations, mais leur expression est longue et non linéaire.

On écrit donc pour chaque degré de liberté :

$$
\dot{\theta}=\omega,\qquad \dot{\omega}=f(\theta,\omega).
$$

Entrées de $f$ :

- $\theta_1,\theta_2$ — positions angulaires des deux bras.
- $\omega_1,\omega_2$ — vitesses angulaires correspondantes.
- constantes physiques : masses, longueurs, gravité, (éventuellement) coefficients d'amortissement.

Sorties de $f$ :

- $\dot{\omega}_1,\dot{\omega}_2$ — accélérations angulaires.

Pourquoi $f$ est-elle compliquée ?

- Les deux pendules sont couplés : le mouvement de l'un modifie les forces et les couples appliqués à l'autre. Cela engendre des termes trigonométriques (sinus, cosinus), des produits de vitesses et des dépendances croisées entre variables.

Comment la simulation numérique l'utilise :

- À chaque pas de temps, l'algorithme évalue $f$ pour calculer les accélérations, puis met à jour les vitesses ($\omega$) et les angles ($\theta$). Cette boucle — évaluer $f$, avancer l'état — est répétée et reconstruit le mouvement du double pendule.

Remarque d'implémentation : on place ces équations sous la forme d'un vecteur d'état

$$
y=[\theta_1,\omega_1,\theta_2,\omega_2]^\top
$$

et on passe à l'intégrateur la fonction

$$
F(t,y)=\begin{bmatrix}\omega_1\\[4pt] f_1(\theta_1,\theta_2,\omega_1,\omega_2)\\[4pt] \omega_2\\[4pt] f_2(\theta_1,\theta_2,\omega_1,\omega_2)\end{bmatrix}
$$

Les solveurs Runge–Kutta (p.ex. `Tsit5`) évaluent plusieurs fois $F$ par pas pour obtenir une mise à jour précise.

### Méthodes de Runge–Kutta (RK) — fonctionnement et avantages

Les méthodes de Runge–Kutta forment une famille d'intégrateurs à un pas qui évaluent la dérivée droite $f(t,y)$ plusieurs fois à l'intérieur d'un pas pour construire une mise à jour de plus haut ordre. En comparaison :

- Euler explicite : $y_{n+1}=y_n + h f(t_n,y_n)$ — méthode d'ordre 1 (erreur globale $O(h)$).
- RK classiques (p.ex. RK4) : évaluations multiples par pas, ordre 4 (erreur globale $O(h^4)$), bien meilleur rapport précision/coût que Euler pour pas raisonnables.

Points clés :

- Ordre et précision : un schéma RK d'ordre $p$ a une erreur locale $O(h^{p+1})$ et une erreur globale $O(h^p)$. Augmenter l'ordre réduit fortement l'erreur à pas identique.
- Contrôle d'erreur adaptatif : les paires embarquées (Runge–Kutta–Fehlberg, Dormand–Prince, Tsitouras) fournissent une estimation d'erreur par pas et permettent d'ajuster automatiquement $h$ pour respecter une tolérance numérique.
- Stabilité et raideur : les RK explicites sont efficaces pour problèmes non raides (comme le pendule ici), mais pour problèmes raides il faut des méthodes implicites (BDF, implicit RK) pour la stabilité.
- Coût par pas : RK de haut ordre effectue plusieurs évaluations de $f$ par pas (coût plus élevé que Euler par pas) mais autorise des pas beaucoup plus grands pour une précision équivalente, donc coût total souvent inférieur.

Pourquoi préférer RK à Euler pour le pendule double :

- Meilleure précision et stabilité locale, réduisant l'erreur de trajectoire sur les fenêtres courtes nécessaires au calibrage.
- Les schémas à pas adaptatif (p.ex. `Tsit5` utilisé dans `DifferentialEquations.jl`) ajustent $h$ automatiquement et gèrent efficacement les régions où la solution varie rapidement.
- Pour des simulations longues où la conservation d'énergie est cruciale, considérer en plus des intégrateurs symplectiques (p.ex. Stormer–Verlet) ; toutefois, pour l'estimation de paramètres et la prédiction à court terme, les RK adaptatifs offrent un excellent compromis précision/performance.

En pratique, `Tsit5` (un schéma Runge–Kutta explicite d'ordre 5/4 adapté) est souvent recommandé pour ce type de système non raide : il combine précision élevée, estimation d'erreur embarquée et efficacité numérique.

- Une méthode simple comme Euler explicite est insuffisante ici, car les erreurs numériques seraient rapidement amplifiées par le caractère chaotique du système.

## 8. Stratégie d'estimation des paramètres (guidage pratique)

- Étape 1 — Géométrie & angles initiaux : estimer les longueurs $L_i$ à partir de calibrations ou de mesures d'image (le projet utilise déjà des longueurs mesurées). Convertir les coordonnées suivies en mètres avant de convertir en angles.
- Étape 2 — Vitesses angulaires initiales : estimer $\dot\theta$ par différences finies sur les premières images (différence centrée avec filtrage du bruit aide).
- Étape 3 — Masses & vitesses initiales : optimiser les masses (et éventuellement de petites corrections aux vitesses angulaires initiales) sur une courte fenêtre initiale où l'amortissement est minimal (le projet utilise des fenêtres 0–0.5 s ou 0–1 s pour éviter un désaccord chaotique trop prononcé).
- Étape 4 — Amortissement : estimer les coefficients d'amortissement en ajustant sur des fenêtres ultérieures où une décroissance d'énergie est visible (ou en faisant correspondre la perte d'énergie entre images). Essayer d'abord un amortissement visqueux simple.
- Utiliser des fonctions objectif basées sur l'erreur quadratique des angles ou sur l'erreur en position cartésienne (convertir les angles simulés en positions pour une comparaison géométrique directe).

## 9. Statistiques et métriques d'erreur

Ce chapitre rassemble les métriques statistiques utilisées pour évaluer la qualité d'ajustement et la divergence entre les mesures et la simulation.

- Coefficient de détermination $R^2$ : mesure la proportion de variance expliquée par le modèle. Pour une série de mesures $y_i$ et prédictions  $\hat y_i$ :

  $$
  R^2 = 1 - \frac{\sum_i (y_i-\hat y_i)^2}{\sum_i (y_i-\bar y)^2},
  $$

  où $\bar y$ est la moyenne des mesures. $R^2\approx1$ signifie un bon ajustement.
- RMSE (Root Mean Squared Error) : erreur quadratique moyenne

  $$
  \mathrm{RMSE}=\sqrt{\frac{1}{N}\sum_{i=1}^N (y_i-\hat y_i)^2}.
  $$
- MAE (Mean Absolute Error) : erreur absolue moyenne

  $$
  \mathrm{MAE}=\frac{1}{N}\sum_{i=1}^N |y_i-\hat y_i|.
  $$
- Écart‑type pour l'énergie : pour la série d'énergie observée $E_i$ (ou des résidus $E_i-\hat E_i$), l'écart‑type

  $$
  \sigma_E=\sqrt{\frac{1}{N-1}\sum_{i=1}^N (E_i-\overline{E})^2}
  $$

  quantifie la dispersion ; en présence d'amortissement contrôlé on peut mesurer le taux moyen de décroissance et la variabilité autour de la tendance.
- Exposant de Lyapunov (estimation pratique) : mesurer la sensibilité aux conditions initiales en observant la divergence exponentielle moyenne de deux trajectoires proches. Si $\delta(t)$ est la distance entre trajectoires, on estime l'exposant maximal $\lambda$ par

  $$
  \delta(t) \approx C e^{\lambda t} \quad\Rightarrow\quad \lambda \approx \frac{1}{\Delta t} \left\langle\ln\frac{\delta(t+\Delta t)}{\delta(t)}\right\rangle,
  $$

  évalué sur une fenêtre finie et moyenné dans le temps et sur plusieurs paires d'initialisations légèrement perturbées. Pour un double pendule, $\lambda>0$ signale chaos local.

Calcul sur les divergences mesure vs simulation :

- Erreurs point à point : calculer $\mathrm{RMSE}$ et $\mathrm{MAE}$ sur angles ($\theta$) ou positions cartésiennes converties depuis les angles.
- Erreur d'énergie : suivre $E(t)=T(t)+V(t)$ observée et simulée, mesurer $\Delta E(t)=E_{obs}(t)-E_{sim}(t)$, puis quantifier la moyenne, l'écart‑type $\sigma_{\Delta E}$ et la pente temporelle moyenne (taux de perte d'énergie) via une régression linéaire sur $\Delta E(t)$.
- Divergence temporelle : mesurer la distance $d(t)$ entre la trajectoire observée et la simulation (p.ex. norme euclidienne des positions des masses). Tracer $\ln d(t)$ ; une pente positive suggère divergence exponentielle (lié à $\lambda$), tandis qu'une croissance plus lente indique erreur numérique ou modélisation.
- Analyse en fenêtres courtes et ensembles : à cause du caractère sensible du système, calculer les métriques sur courtes fenêtres glissantes et/ou sur un ensemble d'initialisations perturbées pour obtenir des distributions d'erreur (médiane, quartiles) plutôt qu'une seule valeur scalaire.
- Normalisation : pour comparer séries de différentes amplitudes, normaliser les erreurs par la plage ou l'écart‑type des données (p.ex. RMSE normalisé = RMSE / std(y)).

Ces métriques combinées donnent un portrait statistique utile : ajustement moyen (R²), amplitude d'erreur (RMSE/MAE), robustesse énergétique (écart‑type et taux de perte d'énergie), et sensibilité dynamique (exposant de Lyapunov, divergence temporelle). Utiliser des visualisations (erreurs temporelles, histogrammes d'erreur, $\ln d(t)$) pour diagnostiquer la nature des écarts.

## 10. Chaos, validation et incertitude

- Le pendule double montre une sensibilité aux conditions initiales : valider les résultats sur plusieurs fenêtres courtes et rapporter des métriques de compétence (erreur moyenne, RMSE, ou distribution sur un ensemble d'initialisations bruitées).
- Lors de prédictions au‑delà de la fenêtre observée, fournir des bandes de confiance (prévision en ensemble) et éviter une interprétation excessive d'une trajectoire unique au‑delà du temps de Lyapunov.

## 11. Comparaison avec les équations balistiques

- La balistique (mouvement d'un projectile) modélise typiquement la translation du centre de masse d'un corps seul. Les équations canoniques (sans traînée) sont

$$
x''(t)=0,\qquad y''(t)=-g,
$$

qui sont découplées et admettent des solutions analytiques. Ajouter la traînée couple les composantes de vitesse via $v=\|\mathbf v\|$ et nécessite une intégration numérique, mais la dimension d'état reste petite et les équations sont structurellement plus simples que pour les pendules.

Les pendules utilisent des coordonnées généralisées angulaires et produisent des EDOs non linéaires, géométriquement couplées, avec des termes en $\sin(\cdot)$ et des produits de vitesses angulaires. L'espace de phase du pendule double est de dimension supérieure et peut être chaotique ; la balistique est typiquement non chaotique.

## 12. Conseils pratiques et extensions

Filtrer le bruit de suivi avant de différencier pour estimer les vitesses.

Lors de l'optimisation des masses, garder des bornes réalistes (masses positives et dans une plage plausible) pour éviter des ajustements pathologiques.

Si une prédiction à long terme est nécessaire, envisager d'estimer un terme de forcing stochastique ou d'utiliser des prévisions en ensemble pour représenter l'erreur de modèle et l'incertitude des conditions initiales.

Pour la visualisation, comparer à la fois les séries temporelles d'angles et les traces cartésiennes des positions des masses (convertir les angles avec la cinématique de la §1).

---

## Annexe : dérivation symbolique Euler–Lagrange

Cette annexe développe l'énergie cinétique et potentielle, calcule les dérivées partielles nécessaires, et montre la manipulation d'Euler–Lagrange qui conduit aux formules explicites d'accélération angulaire de la §3.

1) Récapitulatif de la cinématique (pour mémoire) :

$$
\mathbf{r}_1=(L_1\sin\theta_1,\;L_1\cos\theta_1),\qquad
\mathbf{r}_2=(L_1\sin\theta_1+L_2\sin\theta_2,\;L_1\cos\theta_1+L_2\cos\theta_2).
$$

2) Vitesses au carré (identités utiles) :

$$
\lVert\dot{\mathbf{r}}_1\rVert^2=L_1^2\dot\theta_1^2,
\lVert\dot{\mathbf{r}}_2\rVert^2=L_1^2\dot\theta_1^2+L_2^2\dot\theta_2^2+2L_1L_2\dot\theta_1\dot\theta_2\cos\delta.
$$

3) Énergie cinétique et potentielle (développées) :

$$
T=\tfrac12 m_1 L_1^2\dot\theta_1^2 + \tfrac12 m_2\big(L_1^2\dot\theta_1^2 + L_2^2\dot\theta_2^2 + 2L_1L_2\dot\theta_1\dot\theta_2\cos\delta\big),
V=(m_1+m_2)gL_1\cos\theta_1 + m_2 g L_2\cos\theta_2,
$$

avec $\delta=\theta_2-\theta_1$.

Combiner les termes de $T$ pour écrire le Lagrangien $\mathcal{L}=T-V$ comme :

$$
\mathcal{L}=\tfrac12(m_1+m_2)L_1^2\dot\theta_1^2 + \tfrac12 m_2 L_2^2\dot\theta_2^2 + m_2 L_1 L_2\dot\theta_1\dot\theta_2\cos\delta - V.
$$

4) Dérivées partielles par rapport aux vitesses :

$$
\frac{\partial\mathcal{L}}{\partial\dot\theta_1}=(m_1+m_2)L_1^2\dot\theta_1 + m_2 L_1 L_2\dot\theta_2\cos\delta,\\
\frac{\partial\mathcal{L}}{\partial\dot\theta_2}=m_2 L_2^2\dot\theta_2 + m_2 L_1 L_2\dot\theta_1\cos\delta.
$$

5) Dérivées temporelles des expressions ci‑dessus (règle du produit et $\dot\delta=\dot\theta_2-\dot\theta_1$) :

$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_1}\right)=(m_1+m_2)L_1^2\ddot\theta_1 + m_2 L_1 L_2\big(\ddot\theta_2\cos\delta - \dot\theta_2\sin\delta\,(\dot\theta_2-\dot\theta_1)\big),\\
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_2}\right)=m_2 L_2^2\ddot\theta_2 + m_2 L_1 L_2\big(\ddot\theta_1\cos\delta - \dot\theta_1\sin\delta\,(\dot\theta_1-\dot\theta_2)\big).
$$

6) Dérivées partielles par rapport aux angles (prendre en compte les signes provenant de $V$ et la dépendance de $\cos\delta$ aux angles) :

$$
\frac{\partial\mathcal{L}}{\partial\theta_1}=m_2 L_1 L_2\dot\theta_1\dot\theta_2\sin\delta + (m_1+m_2)gL_1\sin\theta_1,\\
\frac{\partial\mathcal{L}}{\partial\theta_2}=-m_2 L_1 L_2\dot\theta_1\dot\theta_2\sin\delta + m_2 g L_2\sin\theta_2.
$$

7) Les équations d'Euler–Lagrange $\frac{d}{dt}(\partial_{\dot\theta_i}\mathcal{L})-\partial_{\theta_i}\mathcal{L}=0$ donnent deux équations scalaires. En substituant les expressions des étapes 5 et 6 et en réarrangeant les termes (regroupant $\ddot\theta_1$, $\ddot\theta_2$, les termes gravitationnels et les termes non linéaires en vitesses) on obtient les formes explicites présentées en §3.

Après simplification algébrique, les deux équations peuvent s'écrire dans la forme compacte utilisée dans l'implémentation ; en introduisant le dénominateur

$$
\mathrm{den}=m_1 + m_2\sin^2\delta
$$

on obtient la paire (reproduite depuis la §3) :

$$
\begin{aligned}
\ddot\theta_1 &= \frac{m_2\sin\delta\big(L_1\dot\theta_1^2\cos\delta + L_2\dot\theta_2^2\big) + m_2 g \sin\theta_2\cos\delta - (m_1 + m_2) g \sin\theta_1}{L_1\,\mathrm{den}},\\
\ddot\theta_2 &= \frac{-m_2 L_2\dot\theta_2^2\sin\delta\cos\delta - (m_1 + m_2) L_1\dot\theta_1^2\sin\delta + (m_1 + m_2) g \sin\theta_1\cos\delta - (m_1 + m_2) g \sin\theta_2}{L_2\,\mathrm{den}}.
\end{aligned}
$$

Remarques :

- Les étapes intermédiaires impliquent de regrouper des termes contenant $\sin\delta$ ou $\cos\delta$ et de factoriser des coefficients communs ; l'annexe montre les dérivées clés de sorte que la réarrangement algébrique est direct à vérifier symboliquement (ou avec un CAS).
