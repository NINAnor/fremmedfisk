DROP TABLE IF EXISTS fremmedfisk.oerekyte;
CREATE TABLE fremmedfisk.oerekyte AS
SELECT *,
    CASE
        WHEN (effekt_scenario_2_med_egenspredning >= 30 OR effekt_scenario_1_med_egenspredning >= 30) AND
		    n_intro_avg_scenario_0_med_egenspredning >= 30 THEN effekt_scenario_1_med_egenspredning_prosent
			ELSE 0.0
		END AS effekt_scenario_1_med_egenspredning_prosent_class,
    CASE
        WHEN (effekt_scenario_2_med_egenspredning >= 30 OR effekt_scenario_1_med_egenspredning >= 30) AND
		    n_intro_avg_scenario_0_med_egenspredning >= 30 THEN effekt_scenario_2_med_egenspredning_prosent
			ELSE 0.0
		END AS effekt_scenario_2_med_egenspredning_prosent_class
FROM (
SELECT a.*, nv.verneform, nt.naturtype_tekst AS naturtype, nt."bmvVerdi_tekst" AS naturtype_verdi,
    n_intro_avg_scenario_0_med_egenspredning - n_intro_avg_scenario_2_med_egenspredning AS effekt_scenario_2_med_egenspredning,
    n_intro_avg_scenario_0_med_egenspredning - n_intro_avg_scenario_1_med_egenspredning AS effekt_scenario_1_med_egenspredning,
    n_intro_scenario_0_uten_egenspredning - n_intro_scenario_2_uten_egenspredning AS effekt_scenario_2_uten_egenspredning,
    n_intro_scenario_0_uten_egenspredning - n_intro_scenario_1_uten_egenspredning AS effekt_scenario_1_uten_egenspredning,
	CASE WHEN n_intro_avg_scenario_0_med_egenspredning = 0 THEN 0.0 ELSE
    (n_intro_avg_scenario_0_med_egenspredning - n_intro_avg_scenario_2_med_egenspredning) / n_intro_avg_scenario_0_med_egenspredning END AS effekt_scenario_2_med_egenspredning_prosent,
	CASE WHEN n_intro_avg_scenario_0_med_egenspredning = 0 THEN 0.0 ELSE
    (n_intro_avg_scenario_0_med_egenspredning - n_intro_avg_scenario_1_med_egenspredning) / n_intro_avg_scenario_0_med_egenspredning END AS effekt_scenario_1_med_egenspredning_prosent,
	CASE WHEN n_intro_scenario_0_uten_egenspredning = 0 THEN 0.0 ELSE
    (n_intro_scenario_0_uten_egenspredning - n_intro_scenario_2_uten_egenspredning) / CAST(n_intro_scenario_0_uten_egenspredning AS numeric) END AS effekt_scenario_2_uten_egenspredning_prosent,
	CASE WHEN n_intro_scenario_0_uten_egenspredning = 0 THEN 0.0 ELSE
    (n_intro_scenario_0_uten_egenspredning - n_intro_scenario_1_uten_egenspredning) / CAST(n_intro_scenario_0_uten_egenspredning AS numeric) END AS effekt_scenario_1_uten_egenspredning_prosent
 FROM (SELECT x.*, m.municipality FROM (
    SELECT
        l.geom, l.waterbody, l.area_km2, t.*, CAST('Ørekyte' AS text) AS species, ft."forekomst_oerret", ft."forekomst_harr", ft."forekomst_gjedde", ft."forekomst_oerkyte", ft."forekomst_solabbor" 
		FROM
        (SELECT
            *
        FROM
             (SELECT "waterBodyID",
                max(n_intro) AS n_inro_max_scenario_2_med_egenspredning,
                avg(n_intro) AS n_intro_avg_scenario_2_med_egenspredning,
                min(n_intro) AS n_intro_min_scenario_2_med_egenspredning,
                max(n_intro) - min(n_intro) AS n_intro_range_scenario_2_med_egenspredning,
                max(p_intro) AS p_intro_max_scenario_2_med_egenspredning,
                avg(p_intro) AS p_intro_avg_scenario_2_med_egenspredning,
                min(p_intro) AS p_intro_min_scenario_2_med_egenspredning,
                max(p_intro) - min(p_intro) AS p_intro_range_scenario_2_med_egenspredning
            FROM
                (SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_2_500_simu_600_deg UNION ALL
                SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_2_500_simu_700_deg UNION ALL
                SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_2_500_simu_800_deg) AS emt_a GROUP BY "waterBodyID") AS s2me LEFT JOIN
			(SELECT "waterBodyID",
                max(n_intro) AS n_inro_max_scenario_1_med_egenspredning,
                avg(n_intro) AS n_intro_avg_scenario_1_med_egenspredning,
                min(n_intro) AS n_intro_min_scenario_1_med_egenspredning,
                max(n_intro) - min(n_intro) AS n_intro_range_scenario_1_med_egenspredning,
                max(p_intro) AS p_intro_max_scenario_1_med_egenspredning,
                avg(p_intro) AS p_intro_avg_scenario_1_med_egenspredning,
                min(p_intro) AS p_intro_min_scenario_1_med_egenspredning,
                max(p_intro) - min(p_intro) AS p_intro_range_scenario_1_med_egenspredning
            FROM
                (SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_1_500_simu_600_deg UNION ALL
                SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_1_500_simu_700_deg UNION ALL
                SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_1_500_simu_800_deg) AS emt_a GROUP BY "waterBodyID") AS emt USING ("waterBodyID") LEFT JOIN
            (SELECT "waterBodyID",
                max(n_intro) AS n_inro_max_scenario_0_med_egenspredning,
                avg(n_intro) AS n_intro_avg_scenario_0_med_egenspredning,
                min(n_intro) AS n_intro_min_scenario_0_med_egenspredning,
                max(n_intro) - min(n_intro) AS n_intro_range_scenario_0_med_egenspredning,
                max(p_intro) AS p_intro_max_scenario_0_med_egenspredning,
                avg(p_intro) AS p_intro_avg_scenario_0_med_egenspredning,
                min(p_intro) AS p_intro_min_scenario_0_med_egenspredning,
                max(p_intro) - min(p_intro) AS p_intro_range_scenario_0_med_egenspredning
            FROM
                (SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_0_500_simu_600_deg UNION ALL
                SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_0_500_simu_700_deg UNION ALL
                SELECT "waterBodyID", CAST(n_intro AS smallint), p_intro FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_0_500_simu_800_deg) AS eut_a
                GROUP BY "waterBodyID") AS eut NATURAL LEFT JOIN
            (SELECT "waterBodyID", CAST(n_intro AS smallint) AS n_intro_scenario_2_uten_egenspredning, p_intro AS p_intro_scenario_2_uten_egenspredning FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_2_500_simu_no_second) AS s2ue USING ("waterBodyID") LEFT JOIN
            (SELECT "waterBodyID", CAST(n_intro AS smallint) AS n_intro_scenario_1_uten_egenspredning, p_intro AS p_intro_scenario_1_uten_egenspredning FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_1_500_simu_no_second) AS imt USING ("waterBodyID") LEFT JOIN
            (SELECT "waterBodyID", CAST(n_intro AS smallint) AS n_intro_scenario_0_uten_egenspredning, p_intro AS p_intro_scenario_0_uten_egenspredning FROM fremmedfisk.sim_phoxinus_phoxinus_scenario_0_500_simu_no_second) AS iut USING ("waterBodyID")) AS t NATURAL LEFT JOIN
            (SELECT id AS "waterBodyID", geom, waterbody, area_km2 FROM nofa.lake) AS l LEFT JOIN
			    (SELECT
        *
    FROM
        (SELECT "waterBodyID", CAST(1 AS smallint) AS "forekomst_oerret" FROM nofa.get_last_occurrence_status(CAST(26165 AS integer), countrycodes => CAST('NO' AS text))) AS a NATURAL FULL JOIN
        (SELECT "waterBodyID", CAST(1 AS smallint) AS "forekomst_harr" FROM nofa.get_last_occurrence_status(CAST(26178 AS integer), countrycodes => CAST('NO' AS text))) AS b NATURAL FULL JOIN
        (SELECT "waterBodyID", CAST(1 AS smallint) AS "forekomst_gjedde" FROM nofa.get_last_occurrence_status(CAST(26181 AS integer), countrycodes => CAST('NO' AS text))) AS c NATURAL FULL JOIN
        (SELECT "waterBodyID", CAST(1 AS smallint) AS "forekomst_oerkyte" FROM nofa.get_last_occurrence_status(CAST(26138 AS integer), countrycodes => CAST('NO' AS text))) AS d NATURAL FULL JOIN
        (SELECT "waterBodyID", CAST(1 AS smallint) AS "forekomst_solabbor" FROM nofa.get_last_occurrence_status(CAST(26436 AS integer), countrycodes => CAST('NO' AS text))) AS e) AS ft USING ("waterBodyID")
		) AS x,
        (SELECT geom, municipality FROM "AdministrativeUnits"."Fenoscandia_Municipality_polygon") AS m WHERE ST_Intersects(m.geom, x.geom)) AS a LEFT JOIN
		"Habitats_biotopes"."Norway_ProtectedAreas" AS nv ON ST_Intersects(a.geom, nv.geom) LEFT JOIN
        "Habitats_biotopes"."Norway_NatureTypes" AS nt ON ST_Intersects(a.geom, nt.geom)
		) AS sum_table
		ORDER BY effekt_scenario_1_med_egenspredning DESC;



