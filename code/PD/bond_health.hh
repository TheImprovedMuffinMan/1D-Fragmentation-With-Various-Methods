#ifndef BOND_HEALTH_HH
#define BOND_HEALTH_HH

#include <string>
#include <vector>

/**
 * @brief Computes the remaining health of each bond.
 *
 * A bond health of 1 corresponds to an intact bond, while a bond health of 0
 * corresponds to a fully broken bond.
 *
 * The force in a bond is
 *
 * @f[
 *     F = \mathrm{health} \, k \, s
 * @f]
 *
 * where @f$k@f$ is the bond stiffness and @f$s@f$ is the bond stretch.
 * The health therefore scales the bond stiffness as damage develops.
 *
 * The health depends on the largest stretch the bond has ever reached,
 * so damage is irreversible and bonds cannot heal.
 *
 * Two bond laws are supported:
 *
 * - `"bilinear"`:
 *   The force increases linearly up to the elastic limit @p s_elastic,
 *   then decreases linearly to zero at @p s_break. Between these two
 *   limits, the bond is partially damaged and undergoes softening.
 *
 * @verbatim
 *          force
 *            |      /\
 *            |     /  \
 *            |    /    \
 *            |   /      \
 *            |__/________\______ stretch
 *              0 s_elastic s_break
 * @endverbatim
 *
 * - `"pmb"`:
 *   The bond remains fully elastic until @p s_break and then breaks
 *   completely. This corresponds to Silling's prototype microelastic
 *   brittle bond model.
 *
 * @param law Bond law to use: `"bilinear"` or `"pmb"`.
 * @param max_stretch Maximum stretch ever reached by each bond.
 * @param s_elastic Elastic stretch limit of each bond.
 * @param s_break Breaking stretch of each bond.
 *
 * @return Health of each bond, bounded between 0 and 1.
 */
std::vector<double> bondHealth(const std::string& law,
                               const std::vector<double>& max_stretch,
                               const std::vector<double>& s_elastic,
                               const std::vector<double>& s_break);

#endif