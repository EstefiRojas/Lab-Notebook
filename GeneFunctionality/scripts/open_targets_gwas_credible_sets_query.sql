-- This query finds the single variant with the highest posteriorProbability within each credible set.
-- It also provides counts of variants within that set and credible sets within the parent study.
-- Finally, it filters the results to keep non-colocalizing loci, and for colocalizing loci,
-- it keeps only the one from the study with the most total variants.

-- CTE 1: Select and filter the initial pool of credible sets that meet our criteria.
-- This forms the base for all subsequent calculations.
WITH filtered_credible_sets AS (
  SELECT
    cs.studyId,
    cs.studyLocusId,
    cs.pValueMantissa,
    cs.pValueExponent,
    cs.locus.list AS variants_in_set,
    -- Get the total number of variants for each credible set
    ARRAY_LENGTH(cs.locus.list) AS num_variants_in_set
  FROM
    `bigquery-public-data.open_targets_platform.credible_set` AS cs
  WHERE
    cs.studyType = "gwas"
    -- Filter for p-value < 5e-8
    AND (
      cs.pValueExponent < -8
      OR (cs.pValueExponent = -8 AND cs.pValueMantissa < 5)
    )
    -- Filter for sets where probabilities sum to 1
    -- AND (
    --  SELECT SUM(l.element.posteriorProbability)
    --  FROM UNNEST(cs.locus.list) AS l
    -- ) = 1
),

-- CTE 2: Calculate study-level statistics.
study_level_counts AS (
  SELECT
    studyId,
    COUNT(studyLocusId) AS total_sets_in_study,
    SUM(num_variants_in_set) AS total_variants_in_study
  FROM
    filtered_credible_sets
  GROUP BY
    studyId
),

-- CTE 3: Unnest the variants and rank them by posteriorProbability within each credible set.
ranked_variants AS (
  SELECT
    fcs.studyId,
    fcs.studyLocusId,
    fcs.num_variants_in_set,
    variant.element.variantId,
    fcs.pValueMantissa as pValMant_lead,
    fcs.pValueExponent as pValExp_lead,
    variant.element.pValueMantissa,
    variant.element.pValueExponent,
    variant.element.posteriorProbability,
    ROW_NUMBER() OVER(PARTITION BY fcs.studyLocusId ORDER BY variant.element.posteriorProbability DESC) as rank_num
  FROM
    filtered_credible_sets AS fcs,
    UNNEST(fcs.variants_in_set) AS variant
),

-- CTE 4: Combine the data to get one top variant per credible set, with all study stats.
preliminary_results AS (
  SELECT
    rv.studyId,
    rv.studyLocusId,
    slc.total_sets_in_study,
    slc.total_variants_in_study,
    rv.num_variants_in_set,
    rv.variantId,
    rv.pValMant_lead,
    rv.pValExp_lead,
    rv.pValueMantissa,
    rv.pValueExponent,
    rv.posteriorProbability
  FROM
    ranked_variants AS rv
    JOIN study_level_counts AS slc ON rv.studyId = slc.studyId
  WHERE
    rv.rank_num = 1
),

-- CTE 5: De-duplicate variants that appear in multiple studies.
final_results AS (
  SELECT
    * EXCEPT(final_rank_num)
  FROM (
    SELECT
      *,
      ROW_NUMBER() OVER(PARTITION BY variantId ORDER BY total_variants_in_study DESC) as final_rank_num
    FROM
      preliminary_results
  )
  WHERE
    final_rank_num = 1
),

-- NEW CTE 6: Find the maximum study size for any colocalizing partner of a given locus.
coloc_partner_stats AS (
  SELECT
    fr.studyLocusId,
    -- Find the maximum 'total_variants_in_study' among all direct partners
    MAX(partner_fr.total_variants_in_study) as max_partner_study_variants
  FROM
    final_results AS fr
    -- Join to find partners (locus can be on left or right side)
    JOIN (
      SELECT studyLocusIdLeft as locus, studyLocusIdRight as partner FROM `bigquery-public-data.open_targets_platform.colocalization_coloc`
      UNION ALL
      SELECT studyLocusIdRight as locus, studyLocusIdLeft as partner FROM `bigquery-public-data.open_targets_platform.colocalization_coloc`
    ) AS partners ON fr.studyLocusId = partners.locus
    -- Join back to our results to get the stats for the partner
    JOIN final_results AS partner_fr ON partners.partner = partner_fr.studyLocusId
  GROUP BY
    fr.studyLocusId
)

-- Final SELECT: Join our results with the partner stats and apply the final filter.
SELECT
  fr.*
FROM
  final_results AS fr
  -- LEFT JOIN is crucial: it keeps loci that don't have any colocalizing partners.
  LEFT JOIN coloc_partner_stats AS cps ON fr.studyLocusId = cps.studyLocusId
WHERE
  -- Condition 1: Keep the locus if it has no colocalizing partners (the JOIN resulted in NULL).
  cps.studyLocusId IS NULL
  -- Condition 2: Keep the locus if its study's variant count is greater than or equal to
  -- the variant count of its largest colocalizing partner's study.
  OR fr.total_variants_in_study >= cps.max_partner_study_variants
ORDER BY
  fr.studyId,
  fr.studyLocusId;
