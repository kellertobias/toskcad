//! A bounded ordinary-member graph view, not reconstructed geometry.
use super::*;
use std::fmt::Write;
pub struct MemberDocument<'a, 'input> {
    source: &'a File<'input>,
    member: Member,
    text: String,
    work: usize,
    bytes: usize,
}
impl<'a, 'input> MemberDocument<'a, 'input> {
    pub fn source(&self) -> &'a File<'input> {
        self.source
    }
    pub fn member(&self) -> Member {
        self.member
    }
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy)]
pub(super) struct Roots {
    top: u32,
    context: u32,
    shape: u32,
    roots: [u32; 2],
    header: usize,
}
pub(super) fn roots(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    input: &str,
    w: &mut Ledger<'_>,
) -> Result<Roots> {
    // These are readonly lookups of a root and link already structurally proved
    // by representation::read. Cache once, not once per scene member.
    w.step(entities.len().saturating_mul(2))?;
    let (&top, r) = entities
        .iter()
        .find(|(_, r)| matches!(r.as_slice(),[e]if e.name=="ADVANCED_BREP_SHAPE_REPRESENTATION"))
        .expect("File root validated");
    let context = reference(&r[0].args[2]).expect("validated root context reference");
    let shape = entities
        .values()
        .find_map(|r| match r.as_slice() {
            [
                Entity {
                    name: "SHAPE_DEFINITION_REPRESENTATION",
                    args,
                    ..
                },
            ] => match args.as_slice() {
                [Value::Reference(shape), Value::Reference(rep)]
                    if *rep == top
                        && kind(entities, *shape, "PRODUCT_DEFINITION_SHAPE").is_ok() =>
                {
                    Some(*shape)
                }
                _ => None,
            },
            _ => None,
        })
        .expect("File product link validated");
    // The retained Part21 arena has <=1,048,576 entities. Searching at
    // most len+2 positive IDs therefore finds two free u32 IDs; no nominal
    // max-ID assumption or renumbering of original geometry is needed.
    let mut fresh = [0_u32; 2];
    let mut found = 0;
    let mut next = 1_u32;
    while found < 2 {
        w.step(1)?;
        if !entities.contains_key(&next) {
            fresh[found] = next;
            found += 1;
        }
        next += 1;
    }
    let [root, link] = fresh;
    let header = header(input, w)?.len();
    Ok(Roots {
        top,
        context,
        shape,
        roots: [root, link],
        header,
    })
}
impl<'input> File<'input> {
    /// Build a one-member representation view from ORIGINAL entity component
    /// slices and IDs, preserving all numeric/ref/string tokens. Only complex
    /// record wrapping whitespace and the explicit root membership are new.
    /// Generic public readers can verify an ordinary member without reparsing
    /// unrelated bodies. Mapped items refuse: local prototypes are not world.
    pub fn local_member_document<'a>(
        &'a self,
        item: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MemberDocument<'a, 'input>> {
        valid(options)?;
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(self.work)?;
        let member = *self
            .members
            .iter()
            .find(|m| m.item_entity() == item)
            .ok_or(Error::Unsupported)?;
        if member.mapping_link().is_some() {
            return Err(Error::Unsupported);
        }
        let Roots {
            top,
            context,
            shape,
            ..
        } = self.view;
        let set_bytes = (1 + self.entities.len() / 5)
            .saturating_mul(20 * size_of::<usize>() + 11 * size_of::<u32>());
        let pending_bytes = self.entities.len().saturating_mul(size_of::<u32>());
        let fixed = self
            .bytes
            .saturating_add(size_of::<MemberDocument<'_, '_>>());
        if fixed
            .saturating_add(set_bytes)
            .saturating_add(pending_bytes)
            > options.retained_budget_bytes
        {
            return Err(Error::Limit);
        }
        let mut selected = BTreeSet::new();
        let mut pending = Vec::new();
        pending
            .try_reserve_exact(self.entities.len())
            .map_err(|_| Error::Limit)?;
        for id in [member.local_solid_entity(), context, shape] {
            if selected.insert(id) {
                pending.push(id);
            }
        }
        while let Some(id) = pending.pop() {
            w.step(1)?;
            let records = self
                .entities
                .get(&id)
                .expect("root/ref exists before enqueue");
            for e in records {
                walk(&self.entities, &e.args, &mut selected, &mut pending, &mut w)?;
            }
        }
        drop(pending);
        // Product/context metadata that references the original compound root
        // cannot be reduced to this supported ordinary-member view.
        if selected.contains(&top) {
            return Err(Error::Unsupported);
        }
        let [root, link] = self.view.roots;
        let header = &self.input[..self.view.header];
        let definitions = selected
            .iter()
            .map(|id| {
                let records = &self.entities[id];
                records
                    .iter()
                    .map(|e| e.source.expect("actual original component").len())
                    .sum::<usize>()
                    .saturating_add(64 + records.len())
            })
            .sum::<usize>();
        let estimate = header.len().saturating_add(definitions).saturating_add(256);
        if fixed.saturating_add(set_bytes).saturating_add(estimate) > options.retained_budget_bytes
        {
            return Err(Error::Limit);
        }
        let mut text = String::new();
        text.try_reserve_exact(estimate).map_err(|_| Error::Limit)?;
        let bytes = fixed
            .saturating_add(set_bytes)
            .saturating_add(text.capacity())
            .max(
                fixed
                    .saturating_add(set_bytes)
                    .saturating_add(pending_bytes),
            );
        // Exact-reserve capacity is still admitted by the shared ledger, not an
        // assumption about allocator overcapacity.
        let mut memory = crate::step::insert::resources::Memory::new(
            fixed.saturating_add(set_bytes),
            options.retained_budget_bytes,
        )
        .map_err(Error::from)?;
        memory.add(text.capacity()).map_err(Error::from)?;
        w.step(header.len())?;
        text.push_str(header);
        text.push('\n');
        for id in selected {
            let records = &self.entities[&id];
            write!(text, "#{id}=").expect("String write");
            if records.len() > 1 {
                text.push('(');
            }
            for (i, e) in records.iter().enumerate() {
                if i > 0 {
                    text.push(' ');
                }
                let source = e.source.expect("original component retained");
                w.step(source.len() + 1)?;
                text.push_str(source);
            }
            if records.len() > 1 {
                text.push(')');
            }
            text.push_str(";\n");
        }
        w.step(256)?;
        writeln!(
            text,
            "#{root}=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{}),#{context});",
            member.local_solid_entity()
        )
        .expect("String write");
        writeln!(
            text,
            "#{link}=SHAPE_DEFINITION_REPRESENTATION(#{shape},#{root});"
        )
        .expect("String write");
        text.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        w.step(1)?;
        Ok(MemberDocument {
            source: self,
            member,
            text,
            work: w.used,
            bytes,
        })
    }
}
fn header<'a>(input: &'a str, w: &mut Ledger<'_>) -> Result<&'a str> {
    // File construction already proved balanced comments/strings/parentheses
    // and an actual DATA section. Scan only the header, not the geometry body.
    let b = input.as_bytes();
    let mut i = 0;
    let mut depth = 0usize;
    let mut quoted = false;
    loop {
        w.step(1)?;
        if quoted {
            if b[i] == b'\'' {
                if b.get(i + 1) == Some(&b'\'') {
                    i += 2;
                    continue;
                }
                quoted = false;
            }
            i += 1;
            continue;
        }
        if input[i..].starts_with("/*") {
            let n = input[i + 2..].find("*/").expect("validated comment") + 4;
            w.step(n)?;
            i += n;
            continue;
        }
        if b[i] == b'\'' {
            quoted = true;
            i += 1;
            continue;
        }
        if depth == 0
            && input[i..].starts_with("DATA")
            && b.get(i + 4)
                .is_some_and(|v| !v.is_ascii_alphanumeric() && *v != b'_')
        {
            let mut end = i + 4;
            loop {
                if b[end].is_ascii_whitespace() {
                    w.step(1)?;
                    end += 1;
                } else if input[end..].starts_with("/*") {
                    let n = input[end + 2..].find("*/").expect("validated comment") + 4;
                    w.step(n)?;
                    end += n;
                } else {
                    break;
                }
            }
            if b[end] == b';' {
                return Ok(&input[..end + 1]);
            }
        }
        match b[i] {
            b'(' => depth += 1,
            b')' => depth -= 1,
            _ => {}
        }
        i += 1;
    }
}

#[cfg(test)]
mod tests;
