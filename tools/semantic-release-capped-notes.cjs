// GitHub rejects a release body longer than 125,000 characters. A release
// after a long gap lists every commit, so keep whole lines up to a margin
// below that limit and link the complete comparison instead.
const generator = require("@semantic-release/release-notes-generator");

const GITHUB_RELEASE_BODY_LIMIT = 125000;
const NOTES_LIMIT = 120000;

function compareUrl(context) {
  const repository = context.env?.GITHUB_REPOSITORY;
  const previous = context.lastRelease?.gitTag;
  const next = context.nextRelease?.gitTag;
  if (!repository || !next) return null;
  const server = context.env?.GITHUB_SERVER_URL || "https://github.com";
  return previous ? `${server}/${repository}/compare/${previous}...${next}` : `${server}/${repository}/commits/${next}`;
}

function capNotes(notes, context = {}, limit = NOTES_LIMIT) {
  if (typeof notes !== "string" || notes.length <= limit) return notes;
  const url = compareUrl(context);
  const footer = `\n\n---\n\nRelease notes truncated to fit GitHub's release size limit.${url ? ` Full list of changes: ${url}` : ""}\n`;
  const room = limit - footer.length;
  const cut = notes.lastIndexOf("\n", room);
  return `${notes.slice(0, cut > 0 ? cut : room).trimEnd()}${footer}`;
}

async function generateNotes(pluginConfig, context) {
  return capNotes(await generator.generateNotes(pluginConfig, context), context);
}

module.exports = { generateNotes, capNotes, NOTES_LIMIT, GITHUB_RELEASE_BODY_LIMIT };
